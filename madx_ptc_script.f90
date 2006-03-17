module madx_ptc_script_module
!This module enables the user to execute a PTC script 
!that enables additional functionality that is not possible with MAD-X scripting language
  use madx_keywords
  implicit none
  save
!============================================================================================  
!  PUBLIC INTERFACE
  public                                      :: execscript


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
  
  print*, "Exiting execscript"
  
end subroutine execscript

end module madx_ptc_script_module
