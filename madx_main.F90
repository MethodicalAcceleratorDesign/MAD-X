program madx_main
  use run_madx
  implicit none

#ifdef _INTEL_IFORT_SET_RECL
  include 'mad_recl.fi'
  open(unit=6,RECL=mad_record_length)
#endif

  call madxm
  stop
end program madx_main
