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

subroutine w_ptc_export_xml(filename)
  use ptc_export_xml_module
  implicit none
  integer filename(*)
  call ptc_export_xml(filename)
end subroutine w_ptc_export_xml

subroutine w_ptc_move_to_layout()
  use madx_ptc_module
  implicit none
  call ptc_move_to_layout()
end subroutine w_ptc_move_to_layout

subroutine w_ptc_read_errors()
  use madx_ptc_module
  implicit none
  call ptc_read_errors()
end subroutine w_ptc_read_errors

subroutine w_ptc_refresh_k()
  use madx_ptc_module
  implicit none
  call ptc_refresh_k()
end subroutine w_ptc_refresh_k

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

subroutine w_ptc_twiss(tab_name, summary_name)
  use madx_ptc_twiss_module
  implicit none
  integer tab_name(*)
  integer summary_name(*)
  call ptc_twiss(tab_name, summary_name)
end subroutine w_ptc_twiss

subroutine w_ptc_normal()
  use madx_ptc_normal_module
  implicit none
  call ptc_normal()
end subroutine w_ptc_normal

subroutine w_ptc_moments(no)
  use madx_ptc_distrib_module
  implicit none
  integer no
  call ptc_moments(no)
end subroutine w_ptc_moments

subroutine w_ptc_initmoments()
  use madx_ptc_distrib_module
  implicit none
  call initmoments()
end subroutine w_ptc_initmoments

subroutine w_ptc_dumpmaps()
  use madx_ptc_module
  implicit none
  call ptc_dumpmaps()
end subroutine w_ptc_dumpmaps

subroutine w_ptc_track(max_obs)
  use madx_ptc_module
  USE madx_ptc_track_run_module, ONLY: ptc_track_run
  implicit none
  integer max_obs
  call ptc_track_run(max_obs)
end subroutine w_ptc_track

subroutine w_ptc_trackline(max_obs)
  use madx_ptc_trackline_module
  implicit none
  integer max_obs
  call ptc_trackline(max_obs)
end subroutine w_ptc_trackline

subroutine w_ptc_track_everystep(max_obs)
  use madx_ptc_trackline_module
  implicit none
  integer max_obs
  call ptc_track_everystep(max_obs)
end subroutine w_ptc_track_everystep

subroutine w_ptc_setdebuglevel(level)
  use madx_ptc_intstate_module
  implicit none
  integer level
  call ptc_setdebuglevel(level)
end subroutine w_ptc_setdebuglevel

subroutine w_ptc_setmapdumplevel(level)
  use madx_ptc_intstate_module
  implicit none
  integer level
  call ptc_setmapdumplevel(level)
end subroutine w_ptc_setmapdumplevel

subroutine w_ptc_setmadprint(level)
  use madx_ptc_intstate_module
  implicit none
  integer level
  call ptc_setmadprint(level)
end subroutine w_ptc_setmadprint

subroutine w_ptc_setseed(seed)
  use madx_ptc_intstate_module
  implicit none
  integer seed
  call ptc_setseed(seed)
end subroutine w_ptc_setseed

subroutine w_ptc_enforce6D(level)
  use madx_ptc_intstate_module
  implicit none
  integer level
  call setenforce6D(level)
end subroutine w_ptc_enforce6D

subroutine w_ptc_putbeambeam()
  use madx_ptc_module
  implicit none
  call putbeambeam()
end subroutine w_ptc_putbeambeam

subroutine w_ptc_setaccel_method(method)
  use madx_ptc_intstate_module
  implicit none
  integer method
  call ptc_setaccel_method(method)
end subroutine w_ptc_setaccel_method

subroutine w_ptc_setradiation(method)
  use precision_constants
  use madx_ptc_intstate_module
  implicit none
  integer method
  call ptc_setradiation(method)
end subroutine w_ptc_setradiation

subroutine w_ptc_setmodulation(method)
  use precision_constants
  use madx_ptc_intstate_module
  implicit none
  integer method
  call ptc_setmodulation(method)
end subroutine w_ptc_setmodulation

subroutine w_ptc_setexactmis(method)
  use precision_constants
  use madx_ptc_intstate_module
  implicit none
  integer method
  call ptc_setexactmis(method)
end subroutine w_ptc_setexactmis

subroutine w_ptc_settotalpath(method)
  use precision_constants
  use madx_ptc_intstate_module
  implicit none
  integer method
  call ptc_settotalpath(method)
end subroutine w_ptc_settotalpath

subroutine w_ptc_settime(method)
  use precision_constants
  use madx_ptc_intstate_module
  implicit none
  integer method
  call ptc_settime(method)
end subroutine w_ptc_settime

subroutine w_ptc_setnocavity(method)
  use precision_constants
  use madx_ptc_intstate_module
  implicit none
  integer method
  call ptc_setnocavity(method)
end subroutine w_ptc_setnocavity

subroutine w_ptc_setfringe(method)
  use precision_constants
  use madx_ptc_intstate_module
  implicit none
  integer method
  call ptc_setfringe(method)
end subroutine w_ptc_setfringe

subroutine w_ptc_setenvelope(method)
  use precision_constants
  use madx_ptc_intstate_module
  implicit none
  integer method
  call ptc_setenvelope(method)
end subroutine w_ptc_setenvelope

subroutine w_ptc_setspin(method)
  use precision_constants
  use madx_ptc_intstate_module
  implicit none
  integer method
  call ptc_setspin(method)
end subroutine w_ptc_setspin

subroutine w_ptc_setstochastic(method)
  use precision_constants
  use madx_ptc_intstate_module
  implicit none
  integer method
  call ptc_setstochastic(method)
end subroutine w_ptc_setstochastic

subroutine w_ptc_end()
  use precision_constants
  use madx_ptc_module
  implicit none
  call ptc_end()
end subroutine w_ptc_end

subroutine w_ptc_getnfieldcomp(fibreidx,ncomp,nval)
  use madx_ptc_module
  implicit none
  real(kind(1d0)) nval
  integer fibreidx
  integer ncomp
  call ptc_getnfieldcomp(fibreidx,ncomp,nval)
end subroutine w_ptc_getnfieldcomp

subroutine w_ptc_getsfieldcomp(fibreidx,ncomp,nval)
  use madx_ptc_module
  implicit none
  real(kind(1d0)) nval
  integer fibreidx
  integer ncomp
  call ptc_getsfieldcomp(fibreidx,ncomp,nval)
end subroutine w_ptc_getsfieldcomp

subroutine w_ptc_setfieldcomp(fibreidx)
  use madx_ptc_module
  implicit none
  integer fibreidx
  call ptc_setfieldcomp(fibreidx)
end subroutine w_ptc_setfieldcomp

subroutine w_ptc_setknobvalue(fibre)
  use madx_ptc_knobs_module
  implicit none
  integer fibre(*)
  call setknobvalue(fibre)
end subroutine w_ptc_setknobvalue

subroutine w_ptc_refreshtables()
  use madx_ptc_knobs_module
  implicit none
  call filltables()
end subroutine w_ptc_refreshtables

subroutine w_ptc_addknob(fibre)
  use madx_ptc_knobs_module
  implicit none
  integer fibre(*)
  call addknob(fibre)
end subroutine w_ptc_addknob

subroutine w_ptc_addknob_i(paramn)
  use madx_ptc_knobs_module
  implicit none
  integer paramn(*)
  call addknobi(paramn)
end subroutine w_ptc_addknob_i

subroutine w_ptc_addmoment(x,px,y,py,t,dp,tableIA, columnIA, parametric )
  use madx_ptc_distrib_module
  implicit none
  integer               :: x,px,y,py,dp,t
  integer               :: columnIA(*)
  integer               :: tableIA(*)
  integer               :: parametric
  call addmoment(x,px,y,py,t,dp,tableIA, columnIA, parametric )
end subroutine w_ptc_addmoment

subroutine w_ptc_writeparresults(filename)
  use madx_ptc_knobs_module
  implicit none
  integer filename(*)
  call writeparresults(filename)
end subroutine w_ptc_writeparresults

subroutine w_ptc_printframes(filename)
  use pointer_lattice
  implicit none
  integer filename(*)

  if (ASSOCIATED(my_ering) .eqv. .false.) then
     if (ASSOCIATED(m_u) .eqv. .false.) then
          call aafail('w_ptc_printframes:','A PTC command without universe created. Program stops')
          return
     endif
     my_ering => m_u%end
  endif

  call printframes(filename)
end subroutine w_ptc_printframes

subroutine w_ptc_printlayout_rootm(filename)
  use madx_ptc_eplacement_module
  implicit none
  integer filename(*)
  call printlayout_rootm(filename)
end subroutine w_ptc_printlayout_rootm

subroutine w_ptc_eplacement(elementidx,rf)
  use madx_ptc_eplacement_module
  implicit none
  integer elementidx
  integer rf
  call place_element(elementidx,rf)
end subroutine w_ptc_eplacement

subroutine w_ptc_addpush(tabname, colname, polinomial, monomial)
  use madx_ptc_knobs_module
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

subroutine w_ptc_open_gino(scriptname)
  use madx_ptc_script_module
  implicit none
  integer scriptname(*)
  call execginoscript(scriptname)
end subroutine w_ptc_open_gino
