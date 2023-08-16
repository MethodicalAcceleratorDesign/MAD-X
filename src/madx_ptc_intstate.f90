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
  public                            :: getnocharge
  public                            :: getdebug
  public                            :: getenforce6D
  public                            :: setenforce6D
  public                            :: ptc_setdebuglevel
  public                            :: ptc_setmapdumplevel
  public                            :: ptc_setmadprint
  public                            :: ptc_setseed
  public                            :: ptc_setaccel_method
  public                            :: ptc_setexactmis
  public                            :: ptc_setradiation
  public                            :: ptc_setmodulation
  public                            :: ptc_settotalpath
  public                            :: ptc_settime
  public                            :: ptc_setnocavity
  public                            :: ptc_setnocharge
  public                            :: ptc_setspin
  public                            :: ptc_setstochastic
  public                            :: ptc_setenvelope
  public                            :: ptc_setfringe
  public                            :: printintstate



  private
  !============================================================================================
  !  PRIVATE
  !    data structures

  logical(lp),            public   :: maxaccel  ! switch saying to make the reference particle to fly always on the crest
  logical(lp),            public   :: enforce6D = .false. ! normally 6D is reduced to 4D if no cavities are present
  ! this switch prevents it. It is needed to calcualte  fg R56 in a chicane
  type (internal_state),  private  :: intstate = default0
  integer,                private  :: debug = 1    ! defines debug level
  logical(lp),            private  :: nocharge=.true. ! switch saying to grab the charge from the beam or not.

  !    routines

  !--none--!

  !============================================================================================

contains

  logical(lp) function getmaxaccel()
    implicit none
    getmaxaccel = maxaccel
    return
  end function getmaxaccel
  !____________________________________________________________________________________________

  logical(lp) function getnocharge()
    implicit none
    getnocharge = nocharge
    return
  end function getnocharge
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
    !if (getdebug() > 1)
    intstate = state
    if (associated(c_%no) .and. getdebug() > 1) call print(intstate,6)
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

    !if (getdebug() > 1)
    print *, "Initializing internal state"

    intstate = intst - nocavity0
    call update_states

    if ( associated(c_%no) .and. getdebug() > 1) call print(intstate,6)
  end subroutine initintstate
  !____________________________________________________________________________________________

  subroutine ptc_resetinternalstate
    implicit none

    if (getdebug() > 1) then
        print *, "Setting internal state to DEFAULT0"
    end if

    intstate = default0
    default = intstate
    call update_states

    if (associated(c_%no) .and. getdebug() > 1) call print(intstate,6)

  end subroutine ptc_resetinternalstate
  !____________________________________________________________________________________________


  subroutine ptc_setdebuglevel(level)
    implicit none
    integer     :: level

    if (level > 0) then
        print *, "Setting debug level to", level
    end if
    debug = level

  end subroutine ptc_setdebuglevel
  !____________________________________________________________________________________________


  subroutine ptc_setmapdumplevel(level)
    use s_extend_poly, only : mapdump ! LD:29.03.2019
    implicit none
    integer     :: level

    if (level > 0) then
        print *, "Setting mapdump level to", level
    end if
    mapdump = level

  end subroutine ptc_setmapdumplevel

  !____________________________________________________________________________________________


  subroutine ptc_setmadprint(level)
    use precision_constants, only : madxprint ! LD:13.01.2022
    implicit none
    integer     :: level

    if (level > 0) then
        print *, "Setting madprint level to", level
    end if
    madxprint = level.ne.0

  end subroutine ptc_setmadprint

  !____________________________________________________________________________________________
  subroutine ptc_setseed(seed)
    USE gauss_dis
    implicit none
    integer     :: seed

    if (getdebug() > 0) then
        print *, "Setting seed to", seed
    end if

     CALL gaussian_seed(seed)

  end subroutine ptc_setseed

  !____________________________________________________________________________________________

  subroutine setenforce6D(flag)
    implicit none
    integer     :: flag
    if (flag == 0) then
       if (getdebug() > 1) then
           print *, "Switching off ENFORCE6D"
       end if
       enforce6D = .false.
    else
       if (getdebug() > 1) then
           print *, "Setting ENFORCE6D"
       end if
       enforce6D = .true.
    endif

  end subroutine setenforce6D
  !____________________________________________________________________________________________

  logical(lp) function getenforce6D()
    implicit none

    getenforce6D = enforce6D
!    print *, getenforce6D

  end function getenforce6D

  !____________________________________________________________________________________________

  subroutine ptc_setaccel_method(flag)
    implicit none
    integer     :: flag

    if (flag == 1) then
       if (getdebug() > 1) then
           print *, "Setting MAX ACCEL"
       end if
       maxaccel = .true.
    endif

  end subroutine ptc_setaccel_method
  !____________________________________________________________________________________________


   subroutine ptc_setexactmis(flag)
    implicit none
    integer    :: flag
    !    print *, "Setting the flag"
    !    print *, "And the flag is", flag
    if (flag == 1) then
       if (getdebug() > 1) then
           print *, "Switching ON exact missaligment"
       end if
       always_exactmis=.true.
    else
       if (getdebug() > 1) then
           print *, "Switching OFF exact missaligment"
       end if
       always_exactmis=.false.
    endif
    default = intstate
    call update_states
    if (associated(c_%no) .and. getdebug() > 1) call print(intstate,6)
  end subroutine ptc_setexactmis
  !____________________________________________________________________________________________

  subroutine ptc_setradiation(flag)
    implicit none
    integer    :: flag


    if (flag == 1) then
       if (getdebug() > 1) then
           print *, "Switching ON radiation"
       end if
       intstate = intstate + radiation0
    else
       if (getdebug() > 1) then
           print *, "Switching OFF radiation"
       end if
       intstate = intstate - radiation0
    endif
    default = intstate
    call update_states
    if (associated(c_%no) .and. getdebug() > 1) call print(intstate,6)
  end subroutine ptc_setradiation
  !____________________________________________________________________________________________

  subroutine ptc_setmodulation(flag)
    implicit none
    integer    :: flag


    if (flag == 1) then
       if (getdebug() > 1) then
           print *, "Switching ON modulation"
       end if
       intstate = intstate + modulation0
    else
       if (getdebug() > 1) then
           print *, "Switching OFF modulation"
       end if
       intstate = intstate - modulation0
    endif
    default = intstate
    call update_states
    if (associated(c_%no) .and. getdebug() > 1) call print(intstate,6)
  end subroutine ptc_setmodulation
  !____________________________________________________________________________________________

  subroutine ptc_setspin(flag)
    implicit none
    integer    :: flag

    if (flag == 1) then
       if (getdebug() > 1) then
           print *, "Switching ON spin"
       end if
       intstate = intstate + spin0
    else
       if (getdebug() > 1) then
           print *, "Switching OFF spin"
       end if
       intstate = intstate - spin0
    endif

    default = intstate
    call update_states
    if (associated(c_%no) .and. getdebug() > 1) call print(intstate,6)
  end subroutine ptc_setspin
  !____________________________________________________________________________________________

  subroutine ptc_setstochastic(flag)
    implicit none
    integer    :: flag

    if (flag == 1) then
       if (getdebug() > 1) then
           print *, "Switching ON stochastic"
       end if
       intstate = intstate + stochastic0
    else
       if (getdebug() > 1) then
           print *, "Switching OFF stochastic"
       end if
       intstate = intstate - stochastic0
    endif

    default = intstate
    call update_states
    if (associated(c_%no) .and. getdebug() > 1) call print(intstate,6)
  end subroutine ptc_setstochastic


  !____________________________________________________________________________________________

  subroutine ptc_setenvelope(flag)
    implicit none
    integer    :: flag

    if (flag == 1) then
       if (getdebug() > 1) then
           print *, "Switching ON envelope"
       end if
       intstate = intstate + envelope0
    else
       if (getdebug() > 1) then
           print *, "Switching OFF envelope"
       end if
       intstate = intstate - envelope0
    endif

    default = intstate
    call update_states
    if (associated(c_%no) .and. getdebug() > 1) call print(intstate,6)
  end subroutine ptc_setenvelope

  !____________________________________________________________________________________________

  subroutine ptc_setfringe(flag)
    implicit none
    integer    :: flag

    if (flag == 1) then
       if (getdebug() > 1) then
           print *, "Switching ON fringe"
       end if
       intstate = intstate + fringe0
    else
       if (getdebug() > 1) then
           print *, "Switching OFF fringe"
       end if
       intstate = intstate - fringe0
    endif

    default = intstate
    call update_states
    if (associated(c_%no) .and. getdebug() > 1) call print(intstate,6)
  end subroutine ptc_setfringe

  !____________________________________________________________________________________________

  subroutine ptc_settotalpath(flag)
    implicit none
    integer    :: flag

    if (flag == 1) then
       if (getdebug() > 1) then
           print *, "Switching ON totalpath (and switching OFF delta and only_4d)"
       end if

       !intstate = intstate  + totalpath0
       ! actually this is done automatically by PTC: only_4d can not be together with totalpath, detla can be on exclusively with only_4d
       intstate = intstate - delta0 - only_4d0 + totalpath0
    else
       if (getdebug() > 1) then
           print *, "Switching OFF totalpath"
       end if
       intstate = intstate - totalpath0
    endif

    default = intstate
    call update_states
    if (associated(c_%no) .and. getdebug() > 1) call print(intstate,6)

  end subroutine ptc_settotalpath
  !____________________________________________________________________________________________


  subroutine ptc_settime(flag)
    implicit none
    integer    :: flag

    !    print *, "Setting the flag"
    !    print *, "And the flag is", flag

    if (flag == 1) then
       if (getdebug() > 1) then
           print *, "Switching ON time"
       end if
       intstate = intstate + time0
    else
       if (getdebug() > 1) then
           print *, "Switching OFF time"
       end if
       intstate = intstate - time0
    endif

    default = intstate
    call update_states
    if (associated(c_%no) .and. getdebug() > 1) call print(intstate,6)

  end subroutine ptc_settime

  !____________________________________________________________________________________________

  subroutine ptc_setnocavity(flag)
    implicit none
    integer    :: flag

    if (flag == 1) then
       if (getdebug() > 1) then
           print *, "Switching ON nocavity"
       end if
       intstate = intstate + nocavity0
    else
       if (getdebug() > 1) then
           print *, "Switching OFF nocavity and (also) delta and only_4d"
       end if
       intstate = intstate  - delta0 - only_4d0 - nocavity0
    endif

    default = intstate
    call update_states
    if (associated(c_%no) .and. getdebug() > 1) call print(intstate,6)
  end subroutine ptc_setnocavity

  !____________________________________________________________________________________________

  subroutine ptc_setnocharge(flag)
    implicit none
    integer     :: flag

    if (flag == 1) then
      if (getdebug() > 1) then
          print *, "Switching ON nocharge"
      end if
      nocharge = .true.
    else
      if (getdebug() > 1) then
          print *, "Switching OFF nocharge (using the charge from the beam)"
      end if
      nocharge = .false.
    endif
  end subroutine ptc_setnocharge

  !____________________________________________________________________________________________

  subroutine printintstate(n)
    implicit none
    integer               :: n
    if (associated(c_%no) ) then
      call print(intstate,n)
    else
      write(n,*) 'printintstate: Can not print, PTC is not initialized yet'
    endif
  end subroutine printintstate
  !____________________________________________________________________________________________

end module madx_ptc_intstate_module
