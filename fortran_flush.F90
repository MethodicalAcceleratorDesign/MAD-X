      subroutine call_fortran_flush
#ifdef _NAG
use f90_unix
#endif
      implicit none
      call flush(6)
      return
      end
