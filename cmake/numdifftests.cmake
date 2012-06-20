get_target_property(binaryname madxbin LOCATION)
get_target_property(ndiffbin numdiff LOCATION)

execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink 
   ${CMAKE_SOURCE_DIR}/examples ${CMAKE_BINARY_DIR}/examples)

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
endmacro()

# First parameter is test name,
# second is a bool saying if it is 
# a long test (test-user)

numdiff_test(test-ibs False)
numdiff_test(test-jacobian False)
numdiff_test(test-jacobian-2 False)
numdiff_test(test-jacobian-knobs False)
numdiff_test(test-ptc-twiss False)
numdiff_test(test-ptc-normal False)

# Tests that require afs:
if(EXISTS "/afs/cern.ch/")
   numdiff_test(test-twiss False)
   numdiff_test(test-match True)
else()
   message(STATUS "afs is not available, some tests will be missing")
endif()
