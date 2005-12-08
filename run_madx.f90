module run_madx
  use madx_ptc_module
  implicit none

contains

  subroutine madxm


    !--- Fortran main program to call C control module

    implicit none

    call madx()


  end subroutine madxm

end module run_madx
