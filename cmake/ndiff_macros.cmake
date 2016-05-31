
if(WIN32)
   if(NOT EXISTS ${CMAKE_BINARY_DIR}/tests/share)
      message(STATUS "Copying tests/share folder, this will take some time...")
      execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory
         ${CMAKE_SOURCE_DIR}/tests/share ${CMAKE_BINARY_DIR}/tests/share)
   endif()
else()
   file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests)
   execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink 
      ${CMAKE_SOURCE_DIR}/tests/share ${CMAKE_BINARY_DIR}/tests/share)
endif()

macro(ndiff_test testname islong)
   if(${islong})
      set(_testname ${testname}_LONG)
   else()
      set(_testname ${testname})
   endif()
   execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory
      ${CMAKE_CURRENT_SOURCE_DIR}/tests/${testname}
      ${CMAKE_CURRENT_BINARY_DIR}/tests/${testname})
  add_test(NAME ${_testname}
      COMMAND ${CMAKE_COMMAND}
      -DTEST_PROG=$<TARGET_FILE:madxbin>
      -DSOURCEDIR=${CMAKE_CURRENT_BINARY_DIR}/tests/${testname}
      -DTEST_NAME=${testname}
      -DNDIFF=$<TARGET_FILE:ndiff>
      -P ${CMAKE_SOURCE_DIR}/cmake/ctestbase.cmake)
   set_tests_properties (${_testname}
      PROPERTIES PASS_REGULAR_EXPRESSION ".*${testname}.*PASS")
   if(NOT ${islong})
      # short tests should never be allowed to take longer than 5 seconds!
      set_tests_properties (${_testname}
         PROPERTIES TIMEOUT 10)
   endif()
endmacro()
