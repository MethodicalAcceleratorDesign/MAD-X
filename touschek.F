      subroutine touschek
      implicit none

#ifdef _WIN32
      !DEC$ ATTRIBUTES C, ALIAS:'_ibs_' :: ibs
      include 'win32calls.fi'
#endif

      print*,'Just a dummy'

      end
