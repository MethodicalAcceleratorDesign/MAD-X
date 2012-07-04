get_target_property(binaryname madxbin LOCATION)
get_target_property(ndiffbin numdiff LOCATION)

if(WIN32)
   if(NOT EXISTS ${CMAKE_BINARY_DIR}/examples)
      execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory
         ${CMAKE_SOURCE_DIR}/examples ${CMAKE_BINARY_DIR}/examples)
   endif()
else()
   execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink 
      ${CMAKE_SOURCE_DIR}/examples ${CMAKE_BINARY_DIR}/examples)
endif()

set(BASESCRIPT ${CMAKE_SOURCE_DIR}/cmake/ctestbase.cmake)

macro(numdiff_test testname testout islong)
   if(${islong})
      set(_testname ${testname}_LONG)
   else()
      set(_testname ${testname})
   endif()
   set(_testout "${testname} ${testout}")
   execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory
      ${CMAKE_CURRENT_SOURCE_DIR}/tests/${testname}
      ${CMAKE_CURRENT_BINARY_DIR}/tests/${testname})
   add_test(${_testname}
      ${CMAKE_COMMAND}
      -DTEST_PROG=${binaryname}
      -DTEST_OUTPUT=${_testout}
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

numdiff_test(test-ibs "ibs_output.tfs" 0)
numdiff_test(test-jacobian "" 0)
numdiff_test(test-jacobian-2 "" 0)
numdiff_test(test-jacobian-knobs "knobfile" 0)
numdiff_test(test-ptc-twiss "Maxwellian_bend_for_ptc.txt fort.18 ptc-twiss-table" 0)
numdiff_test(test-ptc-normal "distort_1_f_end distort_1_h_end fc.3 fort.18 ptc_map_table.tfs ptc_normal_results.tfs" 0)
numdiff_test(test-rfmultipole "fc.16 fc.2 fc.34 fc.3 fc.3.aux test1_flat.seq test1_track.obs0001.p0001 test1_track.obs0002.p0001" 0)
numdiff_test(test-rfmultipole-2 "sectormap" 0)
numdiff_test(test-rfmultipole-3 "sectormap" 0)
numdiff_test(test-rfmultipole-4 "sectormap" 0)

# Tests that require afs:
if(EXISTS "/afs/cern.ch/")
   numdiff_test(test-twiss "sample_optics.tfs" 0)
   numdiff_test(test-match "str.ip8.b1.dat twiss.ir8.b1.data" 1)
else()
   message(STATUS "afs is not available, some tests will be missing")
endif()
