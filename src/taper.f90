subroutine taper(orbit0, iter, step, fname, error)
  !----------------------------------------------------------------------*
  ! Purpose:
  !   Calculate tapering values for the current sequence
  ! Called from mad_taper.c
  ! Input:
  !   orbit0(6)   (real)  initial guess
  !   iter        (int)   number of iterations
  !   step        (real)  step size for stepwise tapering
  !   fname       (char)  filename for output 
  ! Output:
  !   error       (int)   error flag (0: OK, else != 0)
  !
  !  Author: Ghislain Roy BE-ABP
  !  Date: 2021-08-25
  !----------------------------------------------------------------------*
  use io_units, only : tapout
  use twiss0fi, only : fundim
  use matrices, only : EYE
  use math_constfi, only : zero
  use taperfi
  use code_constfi
  use name_lenfi
  implicit none

  double precision :: orbit0(6), step, orbitf(6)
  integer :: iter, error
  integer, parameter :: fnlen = 48
  character(len=fnlen) :: fname
  character(13) :: dummyname = "no_taper_file"
  
  integer :: lg
  character(len=name_len) :: sequ_name
  integer, external :: get_string

  double precision :: rt(6,6), tt(6,6,6), opt(fundim)

  RT  = EYE
  TT  = zero
  OPT = zero
    
  iterate = iter
  stepsize = step
  taperflag = .true.

  if (fname(1:13) .ne. dummyname) then
     open(unit=tapout, file=fname)
     inquire(unit=tapout, opened=taperout)
  endif
     
  lg = get_string('sequence ', 'name ', sequ_name)
  
  if (taperout) then
     write(tapout,'(a,a)')           "! Output of TAPER command for sequence ", sequ_name(:lg)
     write(tapout,'(a,i5,a,e20.12)') "! with parameters ITERATE = ",iterate," STEPSIZE =  ", stepsize
     write(tapout,'(a)')             "! "
  endif
  
  if (stepsize .ne. zero) then
     ! in absence of a set rule for closure condition of tapering for closed orbit calculation, use the first turn version
     call tmfrst(orbit0,orbitf,.true.,.true.,RT,TT,error,0,0,0)
  else
     ! for absolute calculation, closure for tapering follows pt closure.
     call tmclor(orbit0,.true.,.true.,1.e-4,OPT,RT,TT,error)
  endif
  
  ! reset parameters for benefit of "TWISS, TAPERING" command
  if (taperout) close(unit=tapout)
  taperout = .false.
  iterate = 3
  stepsize = zero
  taperflag = .false.
  
end subroutine taper

subroutine taperreset(error)
  !----------------------------------------------------------------------*
  ! Purpose:
  !   Reset all tapering values in current sequence
  ! Called from mad_taper.c
  ! Input:
  !   nil
  ! Output:
  !   error       (int)   error flag (0: OK, else != 0)
  !
  !  Author: Ghislain Roy BE-ABP
  !  Date: 2021-08-12
  !----------------------------------------------------------------------*
  use taperfi
  use code_constfi
  implicit none

  integer :: code, node, error
  integer, external :: restart_sequ, advance_node
  double precision, external :: node_value
  
  node = restart_sequ()       !---  start
10 continue                   !---  loop over nodes
  code = node_value('mad8_type ')

  select case (code)
     
  case (code_rbend, code_sbend)
     call store_node_value('k0 ', 0.d0)
     
  case (code_quadrupole, code_sextupole)
     call store_node_value('ktap ', 0.d0)
     
  end select
  
  if (advance_node().ne.0) then
     node=node+1
     goto 10                 ! loop over nodes
  endif

end subroutine taperreset

