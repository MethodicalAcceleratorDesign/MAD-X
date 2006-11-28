module madx_ptc_script_module
  !This module enables the user to execute a PTC script
  !that enables additional functionality that is not possible with MAD-X scripting language
  use madx_keywords
  !  use pointer_lattice
  implicit none
  save
  private

  !============================================================================================
  !  PUBLIC INTERFACE
  public                                      :: execscript, execginoscript


  !============================================================================================
  !  PRIVATE
  !    data structures

  !    routines

  !============================================================================================

contains
  !____________________________________________________________________________________________

  subroutine execscript(scriptname)
    implicit none
    include 'twissa.fi'
    integer   scriptname(*)
    character(48) scriptfilename

    scriptfilename = charconv(scriptname)
    print*, "I am in execsript: Script name is ", scriptfilename
    CALL read_ptc_command77(scriptfilename)
    print*, "Exiting execscript"

  end subroutine execscript

  subroutine execginoscript(scriptname)
    implicit none
    include 'twissa.fi'
    integer   scriptname(*)
    character(48) scriptfilename

    scriptfilename = charconv(scriptname)
    print*, "I am in execginosript: Script name is ", scriptfilename
    CALL gino_ptc_command77(scriptfilename)
    print*, "Exiting execginoscript"

  end subroutine execginoscript

end module madx_ptc_script_module
