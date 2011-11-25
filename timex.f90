subroutine timex(r1)
  implicit none
  real r1,timestart,timenow
  common /mytimes/timestart
  save
  call timest(0.0)
  call cpu_time(timenow)
  r1=timenow-timestart
  return
end subroutine timex
