module madx_ptc_intstate_module
  use madx_keywords
  implicit none
  save
  public

  public                            :: initintstate
  logical, public                   :: maxaccel  ! switch saying to make the reference particle to fly always on the crest
  ! program automatically
  type (internal_state), public     :: intstate

contains

  subroutine initintstate(intst)
    implicit none
    type (internal_state) :: intst
    print *, "Setting internal state"
    call print(intstate,6)
    intstate = intst-nocavity0
    call update_states
    call print(intstate,6)
  end subroutine initintstate


  subroutine ptc_resetinternalstate
    implicit none

    print *, "Setting internal state to DEFAULT0"
    intstate = default0
    call update_states
    !    call print(intstate,6)
  end subroutine ptc_resetinternalstate
  !____________________________________________________________________________________________


  subroutine ptc_setaccel_method(flag)
    implicit none
    integer     :: flag

    !    print *, "Setting the flag"
    !    print *, "And the flag is", flag

    if (flag == 1) then
       print *, "Setting MAX ACCEL"
       maxaccel = .true.
    endif

  end subroutine ptc_setaccel_method
  !____________________________________________________________________________________________


  subroutine ptc_setexactmis(flag)
    implicit none
    logical     :: flag

    !    print *, "Setting the flag"
    !    print *, "And the flag is", flag

    if (flag) then
       print *, "Switching ON exact missaligment"
       intstate = intstate + EXACTMIS0
    else
       print *, "Switching OFF exact missaligment"
       intstate = intstate - EXACTMIS0
    endif
    call update_states
    !    call print(intstate,6)
  end subroutine ptc_setexactmis
  !____________________________________________________________________________________________

  subroutine ptc_setradiation(flag)
    implicit none
    logical     :: flag

    !    print *, "Setting the flag"
    !    print *, "And the flag is", flag

    if (flag) then
       print *, "Switching ON radiation"
       intstate = intstate + radiation0
    else
       print *, "Switching OFF radiation"
       intstate = intstate - radiation0
    endif
    call update_states
    !    call print(intstate,6)
  end subroutine ptc_setradiation
  !____________________________________________________________________________________________

  subroutine ptc_setfringe(flag)
    implicit none
    logical     :: flag

    !    print *, "Setting the flag"
    !    print *, "And the flag is", flag

    if (flag) then
       print *, "Switching ON fringe"
       intstate = intstate + fringe0
    else
       print *, "Switching OFF fringe"
       intstate = intstate - fringe0
    endif
    call update_states
    !    call print(intstate,6)
  end subroutine ptc_setfringe
  !____________________________________________________________________________________________

  subroutine ptc_settotalpath(flag)
    implicit none
    logical     :: flag

    !    print *, "Setting the flag"
    !    print *, "And the flag is", flag

    if (flag) then
       print *, "Switching ON totalpath"
       intstate = intstate + totalpath0
    else
       print *, "Switching OFF totalpath"
       intstate = intstate - totalpath0
    endif
    call update_states
    !    call print(intstate,6)
  end subroutine ptc_settotalpath
  !____________________________________________________________________________________________


  subroutine ptc_settime(flag)
    implicit none
    logical     :: flag

    !    print *, "Setting the flag"
    !    print *, "And the flag is", flag

    if (flag) then
       print *, "Switching ON time"
       intstate = intstate + time0
    else
       print *, "Switching OFF time"
       intstate = intstate - time0
    endif
    call update_states
    !    call print(intstate,6)
  end subroutine ptc_settime


  subroutine ptc_setnocavity(flag)
    implicit none
    logical     :: flag

    !    print *, "Setting the flag"
    !    print *, "And the flag is", flag

    if (flag) then
       print *, "Switching ON nocavity"
       intstate = intstate + nocavity0
    else
       print *, "Switching OFF nocavity"
       intstate = intstate - nocavity0
    endif
    call update_states
    !    call print(intstate,6)
  end subroutine ptc_setnocavity


  !____________________________________________________________________________________________

end module madx_ptc_intstate_module
