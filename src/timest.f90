subroutine timest(r1)
  implicit none
  real r1,timestart
  logical start
  common /mytimes/timestart
  data start /.false./
  save
  if (.not.start) then
     start=.true.
     call cpu_time(timestart)
  endif
  return
end subroutine timest
