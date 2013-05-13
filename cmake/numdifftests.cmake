include(numdiff_macros)

# First parameter is test name.
# Second is optionally additional output file names.
# Third is a bool saying if it is 
# a long test (test-user).
#
# if your test depend on a different test,
# use set_tests_properties() to define it.

#numdiff_test(test-ibs "ibs_output.tfs" 0)
numdiff_test(test-jacobian 0)
numdiff_test(test-jacobian-2 0)
numdiff_test(test-jacobian-knobs 0)

numdiff_test(test-match 0)
numdiff_test(test-match-2 1)

numdiff_test(test-ptc-twiss 1)
numdiff_test(test-ptc-normal 1)
numdiff_test(test-ptc-trackline 0)
numdiff_test(test-ptc-trackline-2 0)

numdiff_test(test-rfmultipole 0)
numdiff_test(test-rfmultipole-2 0)
numdiff_test(test-rfmultipole-3 0)
numdiff_test(test-rfmultipole-4 0)

numdiff_test(test-twiss 1)
numdiff_test(test-twiss-2 0)
numdiff_test(test-twiss-3 0)
numdiff_test(test-twiss-4 1)
numdiff_test(test-twiss-5 1)
numdiff_test(test-twiss-6 0)
numdiff_test(test-twiss-7 0)

numdiff_test(test-aperture 0)

numdiff_test(test-makethin 0)
numdiff_test(test-makethin-2 1)

numdiff_test(test-survey 0)

numdiff_test(test-track 0)
numdiff_test(test-track-2 1)
numdiff_test(test-track-3 0)
numdiff_test(test-track-4 0)
numdiff_test(test-track-5 0)
numdiff_test(test-track-6 0)
numdiff_test(test-track-7 1)
numdiff_test(test-track-8  1)
numdiff_test(test-track-9  0)
numdiff_test(test-track-10 1)
numdiff_test(test-track-11 1)

set_tests_properties(test-track-2_LONG PROPERTIES DEPENDS test-makethin-2_LONG)
set_tests_properties(test-track-7_LONG PROPERTIES DEPENDS test-makethin-2_LONG)
set_tests_properties(test-track-8_LONG PROPERTIES DEPENDS test-makethin-2_LONG)

numdiff_test(test-emit 0)

numdiff_test(test-touschek 1)
numdiff_test(test-touschek-2 1)

numdiff_test(test-setvars_lin 0)
numdiff_test(test-thick-quad 0)

if(USE_GC)
   numdiff_test(test-memory 1)
   set_tests_properties(test-memory_LONG PROPERTIES TIMEOUT 60)
endif()
