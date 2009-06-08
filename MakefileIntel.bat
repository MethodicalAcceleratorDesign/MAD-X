@ECHO OFF
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

REM Visual C++ set up
if exist "C:\Program Files\Microsoft Visual Studio 9.0\VC\bin" set VCROOT="C:\Program Files\Microsoft Visual Studio 9.0\VC\bin"
if exist "C:\Program Files\Microsoft Visual Studio 9.0\VC\lib" set VCLIB="C:\Program Files\Microsoft Visual Studio 9.0\VC\lib"

if not exist %VCROOT%\vcvars32.bat (
echo %VCROOT%\vcvars32.bat not found. 
echo exit
goto :end
exit
) else (
echo Using Visual C++ in %VCROOT%\vcvars32.bat
)

call %VCROOT%\vcvars32.bat

REM goto :link
REM pause

REM /ML option for single-thread NO LONGER AVAILABLE
SET MULTITHREADING=/MT
REM /ML option does not seem to be recognized!!!
REM MULTITHREADING=/MT
cl /c /EHsc %MULTITHREADING% /Tp .\tpsa.cpp


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
cl -c /Zm1000 -D_FULL -D_CATCH_MEM_W -D_WIN32 %WRAP_FLAG% %MULTITHREADING% %MADX%\fortran_wrappers.c
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 -D_INTEL_IFORT_FLUSH %MADX%\fortran_flush.F90
REM in the above, removing /assume:underscore did not help _call_fortran_flush_ and _flush unresolved
)

cl -c /Zm1000 -D_FULL -D_CATCH_MEM_W -D_WIN32 %WRAP_FLAG% %WRAP_FLAG% %MULTITHREADING% %MADX%\madxp.c
cl -c /Zm1000 -D_WIN32 %WRAP_FLAG% %MULTITHREADING% %MADX%\gxx11psc.c
cl -c /Zm1000 -D_WIN32 %WRAP_FLAG% %MULTITHREADING% %MADX%\rplot.c
cl -c /Zm1000 -D_WIN32 %WRAP_FLAG% %MULTITHREADING% %MADX%\matchptcknobs.c

ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %MADX%\plot.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %MADX%\sodd.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %MADX%\gxx11ps.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %MADX%\dynap.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %MADX%\emit.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %MADX%\twiss.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %MADX%\match.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %MADX%\matchjc.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %MADX%\matchlib.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %MADX%\matchsa.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %MADX%\touschek.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %MADX%\poisson.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %MADX%\survey.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %MADX%\trrun.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %MADX%\util.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %MADX%\orbf.f90
REM ifort /c /names:lowercase /assume:underscore /fpp /O2 %MADX%\ptc_dummy.f90
REM the above fails - ignored
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %MADX%\ibsdb.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %MADX%\resindex.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %FPP%\a_scratch_size.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %FPP%\b_da_arrays_all.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %FPP%\c_dabnew_berz.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %FPP%\c_tpsa_interface.F90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %FPP%\d_lielib.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %FPP%\h_definition.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %FPP%\i_tpsa.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %FPP%\j_tpsalie.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %FPP%\k_tpsalie_analysis.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %FPP%\l_complex_taylor.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %FPP%\m_real_polymorph.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %FPP%\n_complex_polymorph.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %FPP%\o_tree_element.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %FPP%\Sa_extend_poly.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %FPP%\Sb_sagan_pol_arbitrary.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %FPP%\Sc_euclidean.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %FPP%\Sd_frame.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %FPP%\Se_status.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %FPP%\Sf_def_all_kinds.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %FPP%\Sg_sagan_wiggler.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %FPP%\Sh_def_kind.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %FPP%\Si_def_element.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %FPP%\Sk_link_list.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %FPP%\Sl_family.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %FPP%\Sm_tracking.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %FPP%\Sma0_beam_beam_ptc.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %FPP%\Sma_multiparticle.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %FPP%\Sn_mad_like.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %FPP%\So_fitting.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %FPP%\Sp_keywords.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %FPP%\Spb_fake_gino_sub.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %FPP%\Sq_orbit_ptc.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %FPP%\Sqb_accel_ptc.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %FPP%\Sr_spin.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %FPP%\Sra_fitting.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %MADX%\madx_ptc_intstate.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %MADX%\madx_ptc_setcavs.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %MADX%\madx_ptc_knobs.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %MADX%\madx_ptc_module.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %FPP%\St_pointers.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %MADX%\madx_ptc_script.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %MADX%\madx_ptc_eplacement.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %MADX%\madx_ptc_distrib.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %MADX%\madx_ptc_normal.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %MADX%\madx_ptc_twiss.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %MADX%\madx_ptc_track_run.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %MADX%\madx_ptc_trackcavs.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %MADX%\ptc_export_xml.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %MADX%\wrap.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %MADX%\user2_photon.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %MADX%\timest.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %MADX%\timex.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %MADX%\run_madx.f90
ifort /c /names:lowercase /assume:underscore /assume:noold_unit_star -D_INTEL_IFORT_SET_RECL /fpp /O2 %MADX%\madx_main.F90
goto :link
:link

ifort -o %user_main%\madx madx_main.obj run_madx.obj madxp.obj matchptcknobs.obj rplot.obj dynap.obj emit.obj twiss.obj match.obj matchjc.obj matchlib.obj matchsa.obj touschek.obj survey.obj trrun.obj util.obj orbf.obj ibsdb.obj resindex.obj plot.obj sodd.obj gxx11ps.obj gxx11psc.obj a_scratch_size.obj b_da_arrays_all.obj c_dabnew_berz.obj c_tpsa_interface.obj d_lielib.obj h_definition.obj i_tpsa.obj j_tpsalie.obj k_tpsalie_analysis.obj l_complex_taylor.obj m_real_polymorph.obj n_complex_polymorph.obj o_tree_element.obj Sa_extend_poly.obj Sb_sagan_pol_arbitrary.obj Sc_euclidean.obj Sd_frame.obj Se_status.obj Sf_def_all_kinds.obj Sg_sagan_wiggler.obj Sh_def_kind.obj Si_def_element.obj Sk_link_list.obj Sl_family.obj Sm_tracking.obj Sma0_beam_beam_ptc.obj Sma_multiparticle.obj Sn_mad_like.obj So_fitting.obj Sp_keywords.obj Spb_fake_gino_sub.obj Sq_orbit_ptc.obj Sqb_accel_ptc.obj Sr_spin.obj Sra_fitting.obj madx_ptc_intstate.obj madx_ptc_setcavs.obj madx_ptc_knobs.obj madx_ptc_module.obj St_pointers.obj madx_ptc_script.obj madx_ptc_eplacement.obj madx_ptc_distrib.obj madx_ptc_normal.obj madx_ptc_twiss.obj madx_ptc_track_run.obj madx_ptc_trackcavs.obj ptc_export_xml.obj wrap.obj user2_photon.obj poisson.obj timest.obj timex.obj tpsa.obj %WRAPPERS_OBJ%

:end
