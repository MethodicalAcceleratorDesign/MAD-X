      subroutine w_ptc_create_universe()
      implicit none
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC.')
      end

      subroutine w_ptc_create_layout()
      implicit none
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC.')
      end

      subroutine w_ptc_move_to_layout()
      implicit none
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC.')
      end

      subroutine w_ptc_input()
      implicit none
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC.')
      end

      subroutine w_ptc_align()
      implicit none
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC.')
      end

      subroutine w_ptc_track()
      implicit none
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC.')
      end

      subroutine w_ptc_twiss()
      implicit none
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC.')
      end

      subroutine w_ptc_normal()
      implicit none
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC.')
      end

      subroutine w_ptc_moments(no)
      implicit none
      integer no
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC.')
      end

      subroutine w_ptc_initmoments()
      implicit none
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC.')
      end

      subroutine w_ptc_trackline()
      implicit none
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC.')
      end

      subroutine w_ptc_track_everystep()
      implicit none
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC.')
      end

      subroutine  w_ptc_twiss_linac()
      implicit none
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC.')
      end

      subroutine w_ptc_setdebuglevel ()
      implicit none
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC.')
      end

      subroutine w_ptc_enforce6D(level)
      implicit none
      integer level
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC.')
      end subroutine w_ptc_enforce6D

      subroutine w_ptc_setaccel_method ()
      implicit none
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC.')
      end

      subroutine  w_ptc_setexactmis()
      implicit none
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC.')
      end

      subroutine  w_ptc_setradiation()
      implicit none
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC.')
      end

      subroutine  w_ptc_setfringe()
      implicit none
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC.')
      end

      subroutine  w_ptc_settotalpath()
      implicit none
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC.')
      end

      subroutine  w_ptc_settime()
      implicit none
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC.')
      end

      subroutine  w_ptc_setnocavity()
      implicit none
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC.')
      end

      subroutine  w_ptc_addpush()
      implicit none
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC.')
      end

      subroutine  w_ptc_dumpmaps()
      implicit none
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC.')
      end

      subroutine  w_ptc_select()
      implicit none
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC.')
      end subroutine  w_ptc_select

      subroutine  w_ptc_script()
      implicit none
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC.')
      end subroutine  w_ptc_script

      subroutine  w_ptc_open_gino()
      implicit none
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC.')
      end subroutine  w_ptc_open_gino

      subroutine  w_ptc_addknob()
      implicit none
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC.')
      end subroutine  w_ptc_addknob

      subroutine w_ptc_addmoment(x,px,y,py,dp,t,                        &
     &tableIA,columnIA,parametric )
      implicit none
      integer               :: x,px,y,py,dp,t
      integer               :: columnIA(*)
      integer               :: tableIA(*)
      integer               :: parametric
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC.')
      end subroutine w_ptc_addmoment

      integer function w_ptc_getnmoments()
      implicit none
!     Implemented in madx_ptc_distrib.f90
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC.')
      w_ptc_getnmoments = -1
      end function w_ptc_getnmoments

      subroutine  w_ptc_getmomentstabcol(n, tabn, coln)
      implicit none
!     Implemented in madx_ptc_distrib.f90
      integer              :: n
      character(20)        :: tabn
      character(17)        :: coln
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC.')
      end subroutine  w_ptc_getmomentstabcol

      subroutine  w_ptc_addknob_i()
      implicit none
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC.')
      end subroutine  w_ptc_addknob_i

      subroutine  w_ptc_setknobvalue()
      implicit none
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC.')
      end subroutine  w_ptc_setknobvalue

      subroutine  w_ptc_refreshtables()
      implicit none
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC.')
      end subroutine  w_ptc_refreshtables

      subroutine w_ptc_getnfieldcomp(fibreidx,ncomp,nval)
      implicit none
      real nval
      integer fibreidx
      integer ncomp
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC.')
      end subroutine w_ptc_getnfieldcomp

      subroutine w_ptc_getsfieldcomp(fibreidx,ncomp,nval)
      implicit none
      real nval
      integer fibreidx
      integer ncomp
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC.')
      end subroutine w_ptc_getsfieldcomp

      subroutine w_ptc_setfieldcomp(fibreidx)
      implicit none
      integer fibreidx
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC.')
      end subroutine w_ptc_setfieldcomp

      subroutine  w_ptc_printlayout_rootm()
      implicit none
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC.')
      end subroutine  w_ptc_printlayout_rootm

      subroutine  w_ptc_printframes()
      implicit none
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC.')
      end subroutine  w_ptc_printframes

      subroutine  w_ptc_writeparresults()
      implicit none
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC.')
      end subroutine  w_ptc_writeparresults

      subroutine w_ptc_eplacement(elementidx,rf)
      implicit none
      integer elementidx
      integer rf
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC.')
      end subroutine w_ptc_eplacement

      subroutine w_ptc_end()
      implicit none
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC.')
      end subroutine w_ptc_end

      subroutine w_ptc_rviewer()
      implicit none
      call fort_warn('PTC command skipped! ',                           &
     &' Please use proper MAD-X version that includes PTC')
      end subroutine w_ptc_rviewer
