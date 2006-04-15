subroutine w_ptc_create_universe()
  use madx_ptc_module
  implicit none
  call ptc_create_universe()
end subroutine w_ptc_create_universe
subroutine w_ptc_create_layout()
  use madx_ptc_module
  implicit none
  call ptc_create_layout()
end subroutine w_ptc_create_layout
subroutine w_ptc_move_to_layout()
  use madx_ptc_module
  implicit none
  call ptc_move_to_layout()
end subroutine w_ptc_move_to_layout
subroutine w_ptc_input()
  use madx_ptc_module
  implicit none
  call ptc_input()
end subroutine w_ptc_input
subroutine w_ptc_align()
  use madx_ptc_module
  implicit none
  call ptc_align()
end subroutine w_ptc_align

subroutine w_ptc_twiss(tab_name)
  use madx_ptc_module
  implicit none
  integer tab_name(*)
  call ptc_twiss(tab_name)
end subroutine w_ptc_twiss

subroutine w_ptc_dumpmaps()
  use madx_ptc_module
  implicit none
  call ptc_dumpmaps()
end subroutine w_ptc_dumpmaps

subroutine w_ptc_normal()
  use madx_ptc_module
  implicit none
  call ptc_normal()
end subroutine w_ptc_normal

subroutine w_ptc_track(max_obs)
  use madx_ptc_module
  USE madx_ptc_track_run_module, ONLY: ptc_track_run
  implicit none
  integer max_obs
  call ptc_track_run(max_obs)
end subroutine w_ptc_track

subroutine w_ptc_twiss_linac(tab_name)
  use madx_ptc_trackline_module
  implicit none
  integer tab_name(*)
   call ptc_twiss_linac(tab_name)
end subroutine w_ptc_twiss_linac

subroutine w_ptc_trackline(max_obs)
  use madx_ptc_trackline_module
  implicit none
  integer max_obs
  call ptc_trackline(max_obs)
end subroutine w_ptc_trackline

subroutine w_ptc_setdebuglevel(level)
  use madx_ptc_intstate_module
  implicit none
  integer level
  call ptc_setdebuglevel(level)
end subroutine w_ptc_setdebuglevel

subroutine w_ptc_setaccel_method(method)
  use madx_ptc_intstate_module
  implicit none
  integer method
  call ptc_setaccel_method(method)
end subroutine w_ptc_setaccel_method


subroutine w_ptc_setexactmis(method)
  use precision_constants
  use madx_ptc_intstate_module
  implicit none
  logical(lp) method
  call ptc_setexactmis(method)
end subroutine w_ptc_setexactmis

subroutine w_ptc_setradiation(method)
  use precision_constants
  use madx_ptc_intstate_module
  implicit none
  logical(lp) method
  call ptc_setradiation(method)
end subroutine w_ptc_setradiation


subroutine w_ptc_settotalpath(method)
  use precision_constants
  use madx_ptc_intstate_module
  implicit none
  logical(lp) method
  call ptc_settotalpath(method)
end subroutine w_ptc_settotalpath


subroutine w_ptc_settime(method)
  use precision_constants
  use madx_ptc_intstate_module
  implicit none
  logical(lp) method
  call ptc_settime(method)
end subroutine w_ptc_settime

subroutine w_ptc_setnocavity(method)
  use precision_constants
  use madx_ptc_intstate_module
  implicit none
  logical(lp) method
  call ptc_setnocavity(method)
end subroutine w_ptc_setnocavity

subroutine w_ptc_setfringe(method)
  use precision_constants
  use madx_ptc_intstate_module
  implicit none
  logical(lp) method
  call ptc_setfringe(method)
end subroutine w_ptc_setfringe

subroutine w_ptc_end()
  use precision_constants
  use madx_ptc_module
  implicit none
  call ptc_end()
end subroutine w_ptc_end



subroutine w_ptc_addpush(tabname, colname, polinomial, monomial)
  use madx_ptc_tablepush_module
  implicit none
  integer tabname(*)
  integer colname(*)
  integer polinomial
  integer monomial(*)
  
  call addpush(tabname, colname, polinomial, monomial)
    
end subroutine w_ptc_addpush




subroutine w_ptc_script(scriptname)
  use madx_ptc_script_module
  implicit none
  integer scriptname(*)
  
  call execscript(scriptname)
    
end subroutine w_ptc_script
