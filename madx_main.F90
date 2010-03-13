program madx_main
  use run_madx
  implicit none

#ifdef _INTEL_IFORT_SET_RECL
  integer mad_record_length
  parameter (mad_record_length=240)
  open(unit=6,RECL=mad_record_length)
#endif

  call madxm
  stop
end program madx_main
