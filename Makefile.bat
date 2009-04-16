@ECHO OFF
REM User Input - you may need to change the following folder locations
REM This is for the case where WinCVS checks out to F:\MAD\MADXLahey\source
REM
SET user_main=.
REM FPP code
SET FPP=.
REM PTC code
SET PTC=.
REM MAD-X proper code
SET MADX=.

REM on Windows, one should set the 'Include' environment variable from
REM My Computer's Properties.
REM set INCLUDE="C:\Program Files\Microsoft Visual Studio\VC98\include";%INCLUDE%

REM Some environment variables for generating wrappers when calling Fortran from C
IF EXIST C:\Perl\bin\perl.exe (
SET WRAP_FORTRAN_CALLS=1
) ELSE (
SET WRAP_FORTRAN_CALLS=0
)
IF %WRAP_FORTRAN_CALLS%==1 (
ECHO Wrap Fortran calls from C to prevent garbled output.
perl wrap_fortran_calls.pl MSDOS
SET WRAP_FLAG=-D_WRAP_FORTRAN_CALLS
SET WRAPPERS_OBJ=fortran_wrappers.obj fortran_flush.obj
) ELSE (
ECHO No wrapper generated. MAD executable's redirected outputs might be garbled up.
SET WRAP_FLAG=
SET WRAPPERS_OBJ=
)

@ECHO ON

IF %WRAP_FORTRAN_CALLS% == 1 (
cl -c /Zm1000 -D_FULL -D_CATCH_MEM_W -D_WIN32 %WRAP_FLAG% %MADX%\fortran_wrappers.c
lf95 -c -o1 -lfe "-Cpp" -tp %MADX%\fortran_flush.F90
)

cl -c /Zm1000 -D_FULL -D_CATCH_MEM_W -D_WIN32 %WRAP_FLAG% %WRAP_FLAG% %MADX%\madxp.c
cl -c /Zm1000 -D_WIN32 %WRAP_FLAG% %MADX%\gxx11psc.c
cl -c /Zm1000 -D_WIN32 %WRAP_FLAG% %MADX%\rplot.c
cl -c /Zm1000 -D_WIN32 %WRAP_FLAG% %MADX%\matchptcknobs.c

lf95 -c -o1 -tp %MADX%\plot.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %MADX%\sodd.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %MADX%\gxx11ps.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %MADX%\dynap.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %MADX%\emit.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %MADX%\twiss.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %MADX%\match.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %MADX%\matchjc.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %MADX%\matchlib.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %MADX%\matchsa.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %MADX%\touschek.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %MADX%\poisson.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %MADX%\survey.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %MADX%\trrun.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %MADX%\util.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %MADX%\orbf.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %MADX%\ptc_dummy.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %MADX%\ibsdb.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %MADX%\resindex.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %FPP%\a_scratch_size.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %FPP%\b_da_arrays_all.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %FPP%\c_dabnew_berz.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %FPP%\c_tpsa_interface.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %FPP%\d_lielib.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %FPP%\h_definition.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %FPP%\i_tpsa.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %FPP%\j_tpsalie.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %FPP%\k_tpsalie_analysis.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %FPP%\l_complex_taylor.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %FPP%\m_real_polymorph.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %FPP%\n_complex_polymorph.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %FPP%\o_tree_element.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %FPP%\Sa_extend_poly.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %FPP%\Sb_sagan_pol_arbitrary.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %FPP%\Sc_euclidean.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %FPP%\Sd_frame.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %FPP%\Se_status.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %FPP%\Sf_def_all_kinds.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %FPP%\Sg_sagan_wiggler.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %FPP%\Sh_def_kind.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %FPP%\Si_def_element.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %FPP%\Sk_link_list.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %FPP%\Sl_family.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %FPP%\Sm_tracking.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %FPP%\Sma0_beam_beam_ptc.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %FPP%\Sma_multiparticle.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %FPP%\Sn_mad_like.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %FPP%\So_fitting.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %FPP%\Sp_keywords.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %FPP%\Spb_fake_gino_sub.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %FPP%\Sq_orbit_ptc.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %FPP%\Sqb_accel_ptc.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %FPP%\Sr_spin.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %FPP%\Sra_fitting.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %MADX%\madx_ptc_intstate.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %MADX%\madx_ptc_setcavs.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %MADX%\madx_ptc_knobs.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %MADX%\madx_ptc_module.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %FPP%\St_pointers.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %MADX%\madx_ptc_script.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %MADX%\madx_ptc_eplacement.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %MADX%\madx_ptc_distrib.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %MADX%\madx_ptc_normal.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %MADX%\madx_ptc_twiss.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %MADX%\madx_ptc_track_run.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %MADX%\madx_ptc_trackcavs.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %MADX%\ptc_export_xml.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %MADX%\wrap.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %MADX%\user2_photon.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %MADX%\timest.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %MADX%\timex.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %MADX%\run_madx.f90 -winconsole -ml msvc
lf95 -c -o1 -tp %MADX%\madx_main.f90 -winconsole -ml msvc
lf95 -out %user_main%\madx madx_main.obj run_madx.obj madxp.obj matchptcknobs.obj rplot.obj dynap.obj emit.obj twiss.obj match.obj matchjc.obj matchlib.obj matchsa.obj touschek.obj survey.obj trrun.obj util.obj orbf.obj ibsdb.obj resindex.obj plot.obj sodd.obj gxx11ps.obj gxx11psc.obj a_scratch_size.obj b_da_arrays_all.obj c_dabnew_berz.obj c_tpsa_interface.obj d_lielib.obj h_definition.obj i_tpsa.obj j_tpsalie.obj k_tpsalie_analysis.obj l_complex_taylor.obj m_real_polymorph.obj n_complex_polymorph.obj o_tree_element.obj Sa_extend_poly.obj Sb_sagan_pol_arbitrary.obj Sc_euclidean.obj Sd_frame.obj Se_status.obj Sf_def_all_kinds.obj Sg_sagan_wiggler.obj Sh_def_kind.obj Si_def_element.obj Sk_link_list.obj Sl_family.obj Sm_tracking.obj Sma0_beam_beam_ptc.obj Sma_multiparticle.obj Sn_mad_like.obj So_fitting.obj Sp_keywords.obj Spb_fake_gino_sub.obj Sq_orbit_ptc.obj Sqb_accel_ptc.obj Sr_spin.obj Sra_fitting.obj madx_ptc_intstate.obj madx_ptc_setcavs.obj madx_ptc_knobs.obj madx_ptc_module.obj St_pointers.obj madx_ptc_script.obj madx_ptc_eplacement.obj madx_ptc_distrib.obj madx_ptc_normal.obj madx_ptc_twiss.obj madx_ptc_track_run.obj madx_ptc_trackcavs.obj ptc_export_xml.obj wrap.obj user2_photon.obj poisson.obj timest.obj timex.obj %MADX%\tpsa.lib %WRAPPERS_OBJ%
