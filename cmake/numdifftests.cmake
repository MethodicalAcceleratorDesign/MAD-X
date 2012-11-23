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

#numdiff_test(test-ibs "ibs_output.tfs" 0)
numdiff_test(test-jacobian "" 0)
numdiff_test(test-jacobian-2 "" 0)
numdiff_test(test-jacobian-knobs "knobfile" 0)

numdiff_test(test-match "" 0)
numdiff_test(test-match-2 "str.ip8.b1.dat twiss.ir8.b1.data" 1)

numdiff_test(test-ptc-twiss "Maxwellian_bend_for_ptc.txt fort.18 ptc-twiss-table" 1)
numdiff_test(test-ptc-normal "fort.18 Maxwellian_bend_for_ptc.txt ptc_map_table.tfs ptc_normal_results.tfs" 0)

numdiff_test(test-rfmultipole "fc.16 fc.2 fc.34 fc.3 fc.3.aux test1_flat.seq test1_track.obs0001.p0001 test1_track.obs0002.p0001" 0)
numdiff_test(test-rfmultipole-2 "sectormap" 0)
numdiff_test(test-rfmultipole-3 "sectormap" 0)
numdiff_test(test-rfmultipole-4 "sectormap" 0)

numdiff_test(test-twiss "sample_optics.tfs" 1)
numdiff_test(test-twiss-2 "my_sect_file test-twiss-2 twiss_fv9" 0)
numdiff_test(test-twiss-3 "twiss" 0)
numdiff_test(test-twiss-4 "twiss.b1.data twiss.b2.data" 1)

numdiff_test(test-aperture "ap.tfs" 0)

numdiff_test(test-makethin "" 1)

numdiff_test(test-track "out_test-track_ap_collimatorone" 0)
numdiff_test(test-track-2 "track.obs0001.p0001" 1)
numdiff_test(test-track-3 "out_done out_fone" 0)
numdiff_test(test-track-4 "out_done" 0)
numdiff_test(test-track-5 "out_rellipseone" 0)
numdiff_test(test-track-6 "mytab.tfs twiss1.tfs twiss2.tfs track.obs0001.p0001 track.obs0002.p0001 track.obs0001.p0002  track.obs0002.p0002" 0)
numdiff_test(test-track-7 "track.obs0001.p0001 track.obs0001.p0002 track.obs0001.p0003 twissprb1.1" 1)
set_tests_properties(test-track-2_LONG PROPERTIES DEPENDS test-makethin_LONG)
set_tests_properties(test-track-7_LONG PROPERTIES DEPENDS test-makethin_LONG)
