set INCLUDE=C:\PROGRA~1\Microsoft Visual Studio\VC98\Include
cl -c /Zm1000 -D_FULL -D_CATCH_MEM_W -D_WIN32 madxp.c
cl -c /Zm1000 -D_WIN32 gxx11psc.c
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix -tp plot.F
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix -tp sodd.F
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix -tp gxx11ps.F
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix -tp madxm.F
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix -tp dynap.F
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix -tp emit.F
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix -tp twiss.F
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix -tp match.F
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix -tp matchsa.F
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix -tp touschek.F
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix -tp poisson.F
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix -tp survey.F
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix -tp trrun.F
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix -tp util.F
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix -tp orbf.F
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix -tp ptc_dummy.F
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix -tp ibsdb.F
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix -tp resindex.F
lf95 -c -o1 -tp a_scratch_size.f90
lf95 -c -o1 -tp b_da_arrays_all.f90
lf95 -c -o1 -tp c_dabnew.f90
lf95 -c -o1 -tp d_lielib.f90
lf95 -c -o1 -tp e_define_newda.f90
lf95 -c -o1 -tp f_newda.f90
lf95 -c -o1 -tp g_newLielib.f90
lf95 -c -o1 -tp h_definition.f90
lf95 -c -o1 -tp i_tpsa.f90
lf95 -c -o1 -tp j_tpsalie.f90
lf95 -c -o1 -tp k_tpsalie_analysis.f90
lf95 -c -o1 -tp l_complex_taylor.f90
lf95 -c -o1 -tp m_real_polymorph.f90
lf95 -c -o1 -tp n_complex_polymorph.f90
lf95 -c -o1 -tp o_tree_element.f90
lf95 -c -o1 -tp Sa_extend_poly.f90
lf95 -c -o1 -tp Sb_1_pol_template.f90
lf95 -c -o1 -tp Sb_2_pol_template.f90
lf95 -c -o1 -tp Sb_sagan_pol_arbitrary.f90
lf95 -c -o1 -tp Sc_euclidean.f90
lf95 -c -o1 -tp Sd_frame.f90
lf95 -c -o1 -tp Se_status.f90
lf95 -c -o1 -tp Sf_def_all_kinds.f90
lf95 -c -o1 -tp Sg_0_fitted.f90
lf95 -c -o1 -tp Sg_1_fitted.f90
lf95 -c -o1 -tp Sg_1_template_my_kind.f90
lf95 -c -o1 -tp Sg_2_template_my_kind.f90
lf95 -c -o1 -tp Sg_sagan_wiggler.f90
lf95 -c -o1 -tp Sh_def_kind.f90
lf95 -c -o1 -tp Si_def_element.f90
lf95 -c -o1 -tp Sj_elements.f90
lf95 -c -o1 -tp Sk_link_list.f90
lf95 -c -o1 -tp Sl_family.f90
lf95 -c -o1 -tp Sm_tracking.f90
lf95 -c -o1 -tp Sn_mad_like.f90
lf95 -c -o1 -tp So_fitting.f90
lf95 -c -o1 -tp Sp_keywords.f90
lf95 -c -o1 -tp madx_ptc_module.f90
lf95 -c -o1 -tp madx_ptc_track_run.f90
lf95 -c -o1 -tp user2_photon.f90
lf95 -c -o1 -tp wrap.f90
lf95 -c -o1 -tp timest.f90
lf95 -c -o1 -tp timex.f90
lf95 -out madx madxm.obj madxp.obj dynap.obj emit.obj twiss.obj match.obj matchsa.obj touschek.obj survey.obj trrun.obj util.obj orbf.obj ibsdb.obj resindex.obj ptc_dummy.obj plot.obj sodd.obj gxx11ps.obj gxx11psc.obj timest.obj timex.obj
lf95 -out madxdev madxm.obj madxp.obj dynap.obj emit.obj twiss.obj match.obj matchsa.obj touschek.obj survey.obj trrun.obj util.obj orbf.obj ibsdb.obj resindex.obj plot.obj sodd.obj gxx11ps.obj gxx11psc.obj a_scratch_size.obj b_da_arrays_all.obj c_dabnew.obj d_lielib.obj e_define_newda.obj f_newda.obj g_newLielib.obj h_definition.obj i_tpsa.obj j_tpsalie.obj k_tpsalie_analysis.obj l_complex_taylor.obj m_real_polymorph.obj n_complex_polymorph.obj o_tree_element.obj Sa_extend_poly.obj Sb_1_pol_template.obj Sb_2_pol_template.obj Sb_sagan_pol_arbitrary.obj Sc_euclidean.obj Sd_frame.obj Se_status.obj Sf_def_all_kinds.obj Sg_0_fitted.obj Sg_1_fitted.obj Sg_1_template_my_kind.obj Sg_2_template_my_kind.obj Sg_sagan_wiggler.obj Sh_def_kind.obj Si_def_element.obj Sj_elements.obj Sk_link_list.obj Sl_family.obj Sm_tracking.obj Sn_mad_like.obj So_fitting.obj Sp_keywords.obj madx_ptc_module.obj madx_ptc_track_run.obj user2_photon.obj poisson.obj  wrap.obj timest.obj timex.obj
cl -c /Zm1000 -D_CATCH_MEM_W -D_WIN32 madxp.c
lf95 -out mpars madxm.obj madxp.obj
