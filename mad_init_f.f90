subroutine mad_init_f
    implicit none

#ifdef _INTEL_IFORT_SET_RECL
!    integer mad_record_length
!    parameter (mad_record_length=240)
!    open(unit=6,RECL=mad_record_length)
! ld: let's make it simple
    open(unit=6,RECL=240)
#endif

end subroutine mad_init_f

