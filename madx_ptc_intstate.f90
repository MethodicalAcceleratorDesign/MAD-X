module madx_ptc_intstate_module
  use madx_keywords
  implicit none
  save

  !============================================================================================
  !  PUBLIC INTERFACE

  public                            :: getintstate
  public                            :: setintstate
  public                            :: initintstate
  public                            :: getmaxaccel
  public                            :: getdebug
  public                            :: ptc_setdebuglevel
  public                            :: ptc_setaccel_method
  public                            :: ptc_setexactmis
  public                            :: ptc_setradiation
  public                            :: ptc_settotalpath
  public                            :: ptc_settime
  public                            :: ptc_setnocavity
  public                            :: ptc_setfringe



  private
  !============================================================================================
  !  PRIVATE
  !    data structures

  logical(lp),                public   :: maxaccel  ! switch saying to make the reference particle to fly always on the crest
  type (internal_state),  private  :: intstate = default0
  integer,                private  :: debug = 0    ! defines debug level

  !    routines

    !--none--!

  !============================================================================================

contains

  logical(lp) function getmaxaccel()
    implicit none
    logical(lp)  :: getintstate
    getmaxaccel = maxaccel
    return
  end function getmaxaccel
  !____________________________________________________________________________________________

  type (internal_state) function getintstate()
    implicit none
    !returns the internal state
    getintstate = intstate
    return
  end function getintstate
  !____________________________________________________________________________________________

  subroutine setintstate(state)
    implicit none
    type (internal_state)  :: state
    !sets the internal state
    if (getdebug() > 1) print *, "setintstate: Setting internal state"
    intstate = state
    if (getdebug() > 2) call print(intstate,6)
  end subroutine setintstate
  
  !____________________________________________________________________________________________
  
  integer function getdebug()
    implicit none
    getdebug = debug
  end function getdebug
  !____________________________________________________________________________________________

  subroutine initintstate(intst)
    implicit none
    type (internal_state) :: intst
    
    if (getdebug() > 1) print *, "Initializing internal state"
    
    intstate = intst-nocavity0
    call update_states
    
    if (getdebug() > 2) call print(intstate,6)
  end subroutine initintstate
  !____________________________________________________________________________________________

  subroutine ptc_resetinternalstate
    implicit none

    if (getdebug() > 1) print *, "Setting internal state to DEFAULT0"

    intstate = default0
    call update_states

    if (getdebug() > 2) call print(intstate,6)

  end subroutine ptc_resetinternalstate
  !____________________________________________________________________________________________


  subroutine ptc_setdebuglevel(level)
    implicit none
    integer     :: level

    if (level > 0) print *, "Setting debug level to", level
    debug = level

  end subroutine 
  !____________________________________________________________________________________________


  subroutine ptc_setaccel_method(flag)
    implicit none
    integer     :: flag

    if (flag == 1) then
       if (getdebug() > 1) print *, "Setting MAX ACCEL"
       maxaccel = .true.
    endif

  end subroutine ptc_setaccel_method
  !____________________________________________________________________________________________


  subroutine ptc_setexactmis(flag)
    implicit none
    logical(lp)    :: flag

    !    print *, "Setting the flag"
    !    print *, "And the flag is", flag

    if (flag) then
       if (getdebug() > 1) print *, "Switching ON exact missaligment"
       intstate = intstate + EXACTMIS0
    else
       if (getdebug() > 1) print *, "Switching OFF exact missaligment"
       intstate = intstate - EXACTMIS0
    endif
    call update_states
    if (getdebug() > 2) call print(intstate,6)
  end subroutine ptc_setexactmis
  !____________________________________________________________________________________________

  subroutine ptc_setradiation(flag)
    implicit none
    logical(lp)    :: flag


    if (flag) then
       if (getdebug() > 1) print *, "Switching ON radiation"
       intstate = intstate + radiation0
    else
       if (getdebug() > 1) print *, "Switching OFF radiation"
       intstate = intstate - radiation0
    endif
    call update_states
    if (getdebug() > 2) call print(intstate,6)
  end subroutine ptc_setradiation
  !____________________________________________________________________________________________

  subroutine ptc_setfringe(flag)
    implicit none
    logical(lp)    :: flag

    if (flag) then
       if (getdebug() > 1) print *, "Switching ON fringe"
       intstate = intstate + fringe0
    else
       if (getdebug() > 1) print *, "Switching OFF fringe"
       intstate = intstate - fringe0
    endif
    call update_states
    if (getdebug() > 2) call print(intstate,6)
  end subroutine ptc_setfringe
  !____________________________________________________________________________________________

  subroutine ptc_settotalpath(flag)
    implicit none
    logical(lp)    :: flag

    if (flag) then
       if (getdebug() > 1) print *, "Switching ON totalpath"
       intstate = intstate + totalpath0
    else
       if (getdebug() > 1) print *, "Switching OFF totalpath"
       intstate = intstate - totalpath0
    endif

    call update_states
    if (getdebug() > 2) call print(intstate,6)

  end subroutine ptc_settotalpath
  !____________________________________________________________________________________________


  subroutine ptc_settime(flag)
    implicit none
    logical(lp)    :: flag

    !    print *, "Setting the flag"
    !    print *, "And the flag is", flag

    if (flag) then
       if (getdebug() > 1) print *, "Switching ON time"
       intstate = intstate + time0
    else
       if (getdebug() > 1) print *, "Switching OFF time"
       intstate = intstate - time0
    endif

    call update_states
    if (getdebug() > 2) call print(intstate,6)

  end subroutine ptc_settime


  subroutine ptc_setnocavity(flag)
    implicit none
    logical(lp)    :: flag

    if (flag) then
       if (getdebug() > 1) print *, "Switching ON nocavity"
       intstate = intstate + nocavity0
    else
       if (getdebug() > 1) print *, "Switching OFF nocavity"
       intstate = intstate - nocavity0
    endif
    call update_states
    if (getdebug() > 2) call print(intstate,6)
  end subroutine ptc_setnocavity


  !____________________________________________________________________________________________

end module madx_ptc_intstate_module
