      subroutine call_fortran_flush(s)
#ifdef _NAG
use f90_unix
#endif
      implicit none
      character(*) s
      if (len(s) .gt. 0) then
            write(6,*) "flush unit 6 for ", s ! for debug purposes
      endif
      call flush(6)
      return
      end
