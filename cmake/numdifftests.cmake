get_target_property(binaryname madxbin LOCATION)
get_target_property(ndiffbin numdiff LOCATION)

if(WIN32)
   if(NOT EXISTS ${CMAKE_BINARY_DIR}/examples)
      message(STATUS "Copying examples folder, this will take some time...")
      execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory
         ${CMAKE_SOURCE_DIR}/examples ${CMAKE_BINARY_DIR}/examples)
   endif()
   if(NOT EXISTS ${CMAKE_BINARY_DIR}/tests/share)
      message(STATUS "Copying tests/share folder, this will take some time...")
      execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory
         ${CMAKE_SOURCE_DIR}/tests/share ${CMAKE_BINARY_DIR}/tests/share)
   endif()
else()
   execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink 
      ${CMAKE_SOURCE_DIR}/examples ${CMAKE_BINARY_DIR}/examples)
   file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests)
   execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink 
      ${CMAKE_SOURCE_DIR}/tests/share ${CMAKE_BINARY_DIR}/tests/share)
endif()

set(BASESCRIPT ${CMAKE_SOURCE_DIR}/cmake/ctestbase.cmake)

macro(numdiff_test testname islong)
   if(${islong})
      set(_testname ${testname}_LONG)
   else()
      set(_testname ${testname})
   endif()
   execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory
      ${CMAKE_CURRENT_SOURCE_DIR}/tests/${testname}
      ${CMAKE_CURRENT_BINARY_DIR}/tests/${testname})
   add_test(${_testname}
      ${CMAKE_COMMAND}
      -DTEST_PROG=${binaryname}
      -DSOURCEDIR=${CMAKE_CURRENT_BINARY_DIR}/tests/${testname}
      -DTEST_NAME=${testname}
      -DNUMDIFF=${ndiffbin}
      -P ${BASESCRIPT})
   set_tests_properties (${_testname}
      PROPERTIES PASS_REGULAR_EXPRESSION ".*${testname}.*PASS")
endmacro()

# First parameter is test name.
# Second is optionally additional output file names.
# Third is a bool saying if it is 
# a long test (test-user).

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

numdiff_test(test-makethin 1)

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

set_tests_properties(test-track-2_LONG PROPERTIES DEPENDS test-makethin_LONG)
set_tests_properties(test-track-7_LONG PROPERTIES DEPENDS test-makethin_LONG)
set_tests_properties(test-track-8_LONG PROPERTIES DEPENDS test-makethin_LONG)

numdiff_test(test-emit 0)

numdiff_test(test-touschek 1)
numdiff_test(test-touschek-2 1)

numdiff_test(test-setvars_lin 0)
numdiff_test(test-thick-quad 0)
