
# list of c source files
set(csrcfiles madxp.c gxx11c.c matchptcknobs.c rplot.c )
# list of fortran source files
set(fsrcfiles madx_main.F90 Sqb_accel_ptc.f90 Sr_spin.f90 Sra_fitting.f90 madx_ptc_intstate.f90 madx_ptc_setcavs.f90 madx_ptc_knobs.f90 madx_ptc_module.f90 madx_ptc_distrib.f90 madx_ptc_eplacement.f90 madx_ptc_normal.f90 madx_ptc_script.f90 madx_ptc_trackcavs.f90 madx_ptc_track_run.f90 madx_ptc_twiss.f90 match.f90 matchjc.f90 matchlib.f90 matchsa.f90 orbf.f90 plot.f90 poisson.f90 ptc_export_xml.f90 resindex.f90 run_madx.f90 sodd.f90 Spb_fake_gino_sub.f90 St_pointers.f90 survey.f90 timest.f90 timex.f90 touschek.f90 trrun.f90 twiss.f90 user2_photon.f90 wrap.f90 fortran_flush.F90 a_scratch_size.f90 b_da_arrays_all.f90 d_lielib.f90 util.f90 dynap.f90 emit.f90 gxx11.f90 h_definition.f90 ibsdb.f90 i_tpsa.f90 j_tpsalie.f90 k_tpsalie_analysis.f90 l_complex_taylor.f90 m_real_polymorph.f90 n_complex_polymorph.f90 o_tree_element.f90 Sa_extend_poly.f90 Sb_sagan_pol_arbitrary.f90 Sc_euclidean.f90 Sd_frame.f90 Se_status.f90 Sf_def_all_kinds.f90 Sg_sagan_wiggler.f90 Sh_def_kind.f90 Si_def_element.f90 Sk_link_list.f90 Sl_family.f90 Sm_tracking.f90 Sma0_beam_beam_ptc.f90 Sma_multiparticle.f90 Sn_mad_like.f90 So_fitting.f90 Sp_keywords.f90 Sq_orbit_ptc.f90)


# add source files according to NTPSA option...
if (MADX_NTPSA )
  message("NTPSA turned on")
 set(fsrcfiles ${fsrcfiles} c_dabnew_berz.f90 c_tpsa_interface.F90)
 set(csrcfiles ${csrcfiles} tpsa.cpp)
else (MADX_NTPSA )
 set(fsrcfiles ${fsrcfiles} c_dabnew.f90)
endif  (MADX_NTPSA )

#execute python wrapper scripts (you need to be dependent on one of the output files or else this command will never be ran):
# Unsure about dependencies.. Might be an overkill this one.
ADD_CUSTOM_COMMAND(
  OUTPUT c_wrappers.c c_wrappers.h c_prototypes.h c_wrappers_prototypes.h
  DEPENDS ${csrcfiles} ${fsrcfiles}
  COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/*.c . && cp ${CMAKE_CURRENT_SOURCE_DIR}/*.f90 . && python ${CMAKE_CURRENT_SOURCE_DIR}/wrap_C_calls.py
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Creating C wrapper files"
  )
ADD_CUSTOM_COMMAND(
  OUTPUT fortran_wrappers.c fortran_wrappers.h fortran_prototypes.h fortran_wrappers_prototypes.h
  DEPENDS ${fsrcfiles} c_wrappers.c
  COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/wrap_fortran_calls.py && cp ${CMAKE_CURRENT_SOURCE_DIR}/*.h .
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Creating fortran wrapper files"
  )
# execute python wrap_fortran_calls.py
# COMMAND python wrap_C_calls.py


# main source files... 
# TODO: please check this list at some point!
set(srcfiles ${csrcfiles} ${fsrcfiles} fortran_wrappers.c c_wrappers.c )