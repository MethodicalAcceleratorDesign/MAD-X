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
subroutine w_ptc_normal()
  use madx_ptc_module
  implicit none
  call ptc_normal()
end subroutine w_ptc_normal
subroutine w_ptc_track()
  use madx_ptc_module
  implicit none
  call ptc_track()
end subroutine w_ptc_track
subroutine w_ptc_end()
  use madx_ptc_module
  implicit none
  call ptc_end()
end subroutine w_ptc_end
