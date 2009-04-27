      subroutine call_fortran_flush(s)
#ifdef _NAG
use f90_unix
#endif
      implicit none
#ifdef _INTEL_IFORT_FLUSH
      logical exists
#endif
      character(*) s
      if (len(s) .gt. 0) then
            write(6,*) "flush unit 6 for ", s ! for debug purposes
      endif
#ifdef _INTEL_IFORT_FLUSH
      inquire(unit=6, exist = exists)
      if (exists) then
        flush(6)
      endif
#else
      call flush(6)
#endif
      return
      end
