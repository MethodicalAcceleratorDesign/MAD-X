get_target_property(binaryname madxbin LOCATION)
get_target_property(ndiffbin numdiff LOCATION)

set(BASESCRIPT ${CMAKE_SOURCE_DIR}/cmake/ctestbase.cmake)
execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink 
   ${CMAKE_SOURCE_DIR}/examples ${CMAKE_CURRENT_BINARY_DIR}/examples)

macro(numdiff_test testname islong)
   if(${islong})
      set(_testname ${testname}_LONG)
   else()
      set(_testname ${testname})
   endif()
   execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory 
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${CMAKE_CURRENT_BINARY_DIR})
   add_test(${_testname}
       ${CMAKE_COMMAND}
       -DTEST_PROG=${binaryname}
       -DSOURCEDIR=${CMAKE_CURRENT_BINARY_DIR}/${testname}
       -DTEST_NAME=${testname}
       -DNUMDIFF=${ndiffbin}
       -P ${BASESCRIPT})
endmacro()
