      program madxm
!--- Fortran main program to call C control module
      implicit none
#ifdef _WIN32
        !DEC$ ATTRIBUTES C :: madx
#endif

      call madx()
      end
