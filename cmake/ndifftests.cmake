include(ndiff_macros)

# First parameter is test name.
# Second is optionally additional output file names.
# Third is a bool saying if it is 
# a long test (test-user).
#
# if your test depend on a different test,
# use set_tests_properties() to define it.

#ndiff_test(test-ibs "ibs_output.tfs" 0)
ndiff_test(test-jacobian 0)
ndiff_test(test-jacobian-2 0)
ndiff_test(test-jacobian-knobs 0)

ndiff_test(test-match 0)
ndiff_test(test-match-2 1)

ndiff_test(test-ptc-normal 1)
ndiff_test(test-ptc-trackline 0)
ndiff_test(test-ptc-trackline-2 0)
ndiff_test(test-ptc-trackline-3 1)
ndiff_test(test-ptc-track 0)
ndiff_test(test-ptc-track-3 1)

ndiff_test(test-rfmultipole 0)
ndiff_test(test-rfmultipole-2 0)
ndiff_test(test-rfmultipole-3 0)
ndiff_test(test-rfmultipole-4 0)
ndiff_test(test-rfmultipole-5 0)
ndiff_test(test-rfmultipole-6 0)
ndiff_test(test-rfmultipole-ptc-1 0)

ndiff_test(test-twiss 1)
ndiff_test(test-twiss-2 0)
ndiff_test(test-twiss-3 0)
ndiff_test(test-twiss-4 1)
ndiff_test(test-twiss-5 1)
ndiff_test(test-twiss-6 0)
ndiff_test(test-twiss-8 0)
ndiff_test(test-twiss-9 1)
ndiff_test(test-twiss-10 0)
ndiff_test(test-twiss-11 0)
ndiff_test(test-twiss-12 0)
ndiff_test(test-twiss-13 0)

ndiff_test(test-aperture 0)
ndiff_test(test-aperture-2 1)
ndiff_test(test-aperture-3 0)

ndiff_test(test-makethin 0)
ndiff_test(test-makethin-2 1)
ndiff_test(test-makethin-3 1)

ndiff_test(test-survey 0)

ndiff_test(test-track 0)
ndiff_test(test-track-2 1)
ndiff_test(test-track-3 0)
ndiff_test(test-track-4 0)
ndiff_test(test-track-5 0)
ndiff_test(test-track-6 0)
ndiff_test(test-track-7 1)
ndiff_test(test-track-8  1)
ndiff_test(test-track-9  0)
ndiff_test(test-track-10 1)
ndiff_test(test-track-11 1)

set_tests_properties(test-track-2_LONG PROPERTIES DEPENDS test-makethin-2_LONG)
set_tests_properties(test-track-7_LONG PROPERTIES DEPENDS test-makethin-2_LONG)
set_tests_properties(test-track-8_LONG PROPERTIES DEPENDS test-makethin-2_LONG)

ndiff_test(test-emit 0)

ndiff_test(test-touschek 1)
ndiff_test(test-touschek-2 1)

ndiff_test(test-thick-quad 0)
ndiff_test(test-thick-quad-2 0)
ndiff_test(test-thick-quad-3 0)

ndiff_test(test-survey-2 0)

ndiff_test(test-cororbit 1)
ndiff_test(test-cororbit-2 1)
ndiff_test(test-cororbit-3 0)
ndiff_test(test-cororbit-4 1)

ndiff_test(test-emit-2 0)
ndiff_test(test-ibs 0)
ndiff_test(test-ibs-2 0)
ndiff_test(test-ibs-3 0)
ndiff_test(test-ibs-4 1)
ndiff_test(test-error 0)
ndiff_test(test-error-2 0)
ndiff_test(test-error-3 1)
ndiff_test(test-dynap 1)
ndiff_test(test-c6t 1)
ndiff_test(test-c6t-4 1)
set_tests_properties(test-c6t_LONG PROPERTIES DEPENDS test-makethin-2_LONG)
ndiff_test(test-c6t-2 0)
ndiff_test(test-c6t-3 0)

ndiff_test(test-match-3 0)
ndiff_test(test-match-4 0)
ndiff_test(test-match-5 0)
ndiff_test(test-match-6 1)
ndiff_test(test-match-7 1)
ndiff_test(test-match-8 0)

ndiff_test(test-ptc-twiss-1 1)
ndiff_test(test-ptc-twiss-2 1)
ndiff_test(test-ptc-twiss-3 1)
ndiff_test(test-ptc-twiss-4 1)
ndiff_test(test-ptc-twiss-5Dt 1)
ndiff_test(test-ptc-twiss-6D 1)
ndiff_test(test-ptc-twiss-old1 1)
ndiff_test(test-ptc-twiss-old2 1)
ndiff_test(test-ptc-twiss-old3 1)
ndiff_test(test-ptc-twiss-old4 1)
ndiff_test(test-ptc-twiss-old5 1)
if(NOT WIN32)
   ndiff_test(test-ptc-twiss-old6 1)
endif()
ndiff_test(test-ptc-twiss-old7 1)
ndiff_test(test-ptc-twiss-5D 1)
ndiff_test(test-ptc-twiss-56D 1)
ndiff_test(test-ptc-twiss-56Dt 1)
ndiff_test(test-ptc-twiss-56Dl 1)
ndiff_test(test-ptc-twiss-56Dtl 1)
ndiff_test(test-ptc-twiss-56Dt-ini_map_man 1)
ndiff_test(test-ptc-twiss-56Dt-ini_mtx_man 0)
ndiff_test(test-ptc-twiss-56Dt-ini_mtx_tbl 1)
ndiff_test(test-ptc-twiss-normal-genfu 0)
ndiff_test(test-ptc-twiss-normal-5D 1)
ndiff_test(test-ptc-twiss-normal-6D 1)
ndiff_test(test-ptc-twiss-maptable 1)
ndiff_test(test-ptc-twiss-accel-56D 0)
ndiff_test(test-ptc-twiss-6D-ALS 1)


ndiff_test(test-sequence 0)
ndiff_test(test-sequence-2 0)
ndiff_test(test-sequence-3 0)
ndiff_test(test-sequence-4 0)
ndiff_test(test-sequence-5 0)
ndiff_test(test-sequence-6 0)

ndiff_test(test-thick-dipole 0)
ndiff_test(test-thick-dipole-2 0)
ndiff_test(test-thick-dipole-3 0)

ndiff_test(test-plot 1)
ndiff_test(test-plot-2 0)

ndiff_test(test-table 0)
ndiff_test(test-table-2 0)

ndiff_test(test-line 0)

ndiff_test(test-memory 1)
# set_tests_properties(test-memory_LONG PROPERTIES TIMEOUT 60)
ndiff_test(test-beam 0)
set_tests_properties(test-beam PROPERTIES WILL_FAIL 1)

ndiff_test(test-match-9 0)
ndiff_test(test-setknob 0)
ndiff_test(test-fillknob 0)

ndiff_test(test-interpolate 0)
