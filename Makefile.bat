REM User Input - you may need to change the following folder locations
REM This is for the case where WinCVS checks out to F:\MAD\MADXLahey\source
REM
SET user_main=F:\MAD\MADXLahey\source\madX
REM FPP code
SET FPP=F:\MAD\MADXLahey\source\madX
REM PTC code
SET PTC=F:\MAD\MADXLahey\source\madX
REM MAD-X proper code
SET MADX=F:\MAD\MADXLahey\source\madX

set INCLUDE=C:\PROGRA~1\Microsoft Visual Studio\VC98\Include
cl -c /Zm1000 -D_FULL -D_CATCH_MEM_W -D_WIN32 %MADX%\madxp.c
cl -c /Zm1000 -D_WIN32 %MADX%\gxx11psc.c
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix -tp %MADX%\plot.F
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix -tp %MADX%\sodd.F
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix -tp %MADX%\gxx11ps.F
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix -tp %MADX%\dynap.F
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix -tp %MADX%\emit.F
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix -tp %MADX%\twiss.F
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix -tp %MADX%\match.F
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix -tp %MADX%\matchjc.F
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix -tp %MADX%\matchsa.F
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix -tp %MADX%\touschek.F
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix -tp %MADX%\poisson.F
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix -tp %MADX%\survey.F
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix -tp %MADX%\trrun.F
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix -tp %MADX%\util.F
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix -tp %MADX%\orbf.F
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix -tp %MADX%\ptc_dummy.F
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix -tp %MADX%\ibsdb.F
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix -tp %MADX%\resindex.F
lf95 -c -o1 -nconcc -lfe "-Cpp" -fix -tp %MADX%\madxm.F
lf95 -c -o1 -tp %FPP%\a_scratch_size.f90
lf95 -c -o1 -tp %FPP%\b_da_arrays_all.f90
lf95 -c -o1 -tp %FPP%\c_dabnew.f90
lf95 -c -o1 -tp %FPP%\d_lielib.f90
lf95 -c -o1 -tp %FPP%\e_define_newda.f90
lf95 -c -o1 -tp %FPP%\f_newda.f90
lf95 -c -o1 -tp %FPP%\g_newLielib.f90
lf95 -c -o1 -tp %FPP%\h_definition.f90
lf95 -c -o1 -tp %FPP%\i_tpsa.f90
lf95 -c -o1 -tp %FPP%\j_tpsalie.f90
lf95 -c -o1 -tp %FPP%\k_tpsalie_analysis.f90
lf95 -c -o1 -tp %FPP%\l_complex_taylor.f90
lf95 -c -o1 -tp %FPP%\m_real_polymorph.f90
lf95 -c -o1 -tp %FPP%\n_complex_polymorph.f90
lf95 -c -o1 -tp %FPP%\o_tree_element.f90
lf95 -c -o1 -tp %FPP%\Sa_extend_poly.f90
lf95 -c -o1 -tp %FPP%\Sb_1_pol_template.f90
lf95 -c -o1 -tp %FPP%\Sb_2_pol_template.f90
lf95 -c -o1 -tp %FPP%\Sb_sagan_pol_arbitrary.f90
lf95 -c -o1 -tp %FPP%\Sc_euclidean.f90
lf95 -c -o1 -tp %FPP%\Sd_frame.f90
lf95 -c -o1 -tp %FPP%\Se_status.f90
lf95 -c -o1 -tp %FPP%\Sf_def_all_kinds.f90
lf95 -c -o1 -tp %FPP%\Sg_1_fitted.f90
lf95 -c -o1 -tp %FPP%\Sg_1_template_my_kind.f90
lf95 -c -o1 -tp %FPP%\Sg_2_template_my_kind.f90
lf95 -c -o1 -tp %FPP%\Sg_sagan_wiggler.f90
lf95 -c -o1 -tp %FPP%\Sh_def_kind.f90
lf95 -c -o1 -tp %FPP%\Si_def_element.f90
lf95 -c -o1 -tp %FPP%\Sj_elements.f90
lf95 -c -o1 -tp %FPP%\Sk_link_list.f90
lf95 -c -o1 -tp %FPP%\Sl_family.f90
lf95 -c -o1 -tp %FPP%\Sm_tracking.f90
lf95 -c -o1 -tp %FPP%\Sn_mad_like.f90
lf95 -c -o1 -tp %FPP%\So_fitting.f90
lf95 -c -o1 -tp %FPP%\Sp_keywords.f90
lf95 -c -o1 -tp %MADX%\madx_ptc_intstate.f90
lf95 -c -o1 -tp %MADX%\madx_ptc_script.f90
lf95 -c -o1 -tp %MADX%\madx_ptc_setcavs.f90
lf95 -c -o1 -tp %MADX%\madx_ptc_tablepush.f90
lf95 -c -o1 -tp %MADX%\madx_ptc_trackcavs.f90
lf95 -c -o1 -tp %MADX%\madx_ptc_module.f90
lf95 -c -o1 -tp %MADX%\madx_ptc_track_run.f90
lf95 -c -o1 -tp %MADX%\user2_photon.f90
lf95 -c -o1 -tp %MADX%\wrap.f90
lf95 -c -o1 -tp %MADX%\timest.f90
lf95 -c -o1 -tp %MADX%\timex.f90
lf95 -c -o1 -tp %MADX%\run_madx.f90
lf95 -c -o1 -tp %MADX%\madx_main.f90
lf95 -out %user_main%\madx madxm.obj madxp.obj dynap.obj emit.obj twiss.obj match.obj matchjc.obj matchsa.obj touschek.obj survey.obj trrun.obj util.obj orbf.obj ibsdb.obj resindex.obj ptc_dummy.obj plot.obj sodd.obj gxx11ps.obj gxx11psc.obj timest.obj timex.obj
lf95 -out %user_main%\madxdev madx_main.obj run_madx.obj madxp.obj dynap.obj emit.obj twiss.obj match.obj matchjc.obj matchsa.obj touschek.obj survey.obj trrun.obj util.obj orbf.obj ibsdb.obj resindex.obj plot.obj sodd.obj gxx11ps.obj gxx11psc.obj a_scratch_size.obj b_da_arrays_all.obj c_dabnew.obj d_lielib.obj e_define_newda.obj f_newda.obj g_newLielib.obj h_definition.obj i_tpsa.obj j_tpsalie.obj k_tpsalie_analysis.obj l_complex_taylor.obj m_real_polymorph.obj n_complex_polymorph.obj o_tree_element.obj Sa_extend_poly.obj Sb_1_pol_template.obj Sb_2_pol_template.obj Sb_sagan_pol_arbitrary.obj Sc_euclidean.obj Sd_frame.obj Se_status.obj Sf_def_all_kinds.obj Sg_1_fitted.obj Sg_1_template_my_kind.obj Sg_2_template_my_kind.obj Sg_sagan_wiggler.obj Sh_def_kind.obj Si_def_element.obj Sj_elements.obj Sk_link_list.obj Sl_family.obj Sm_tracking.obj Sn_mad_like.obj So_fitting.obj Sp_keywords.obj madx_ptc_module.obj madx_ptc_track_run.obj madx_ptc_intstate.obj madx_ptc_script.obj madx_ptc_setcavs.obj madx_ptc_tablepush.obj madx_ptc_trackcavs.obj user2_photon.obj poisson.obj  wrap.obj timest.obj timex.obj
cl -c /Zm1000 -D_CATCH_MEM_W -D_WIN32 %MADX%\madxp.c
cl -c /Zm1000 -D_CATCH_MEM_W -D_WIN32 %MADX%\rplot.c
lf95 -out %user_main%\mpars madxm.obj madxp.obj
