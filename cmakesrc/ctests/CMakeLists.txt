

if (APPLE)
    set(binaryname "${CMAKE_BINARY_DIR}/${binaryname}.app/Contents/MacOS/./${binaryname}")
else(APPLE)
    set(binaryname "${CMAKE_BINARY_DIR}/./${binaryname}")
endif (APPLE)

ADD_TEST(testruns "testruns.sh" "${binaryname}")
set_tests_properties (testruns 
  PROPERTIES PASS_REGULAR_EXPRESSION "Number of warnings: 0")

# This is a slow test so I commented it out..
# ADD_TEST(testsample "${CMAKE_CURRENT_SOURCE_DIR}/testsample.sh" "${binaryname}" "${CMAKE_CURRENT_SOURCE_DIR}")
# set_tests_properties (testsample 
#   PROPERTIES FAIL_REGULAR_EXPRESSION "Number of errors: [1-9]") # The no warning check fails on this test..

# This is a general way to test easily using a macro script I have written:
# TEST_SCRIPT: name of madx script to execute
# TEST_OUTPUT: name of file madx produces (skips if empty)
add_test(NAME twiss_1
    COMMAND ${CMAKE_COMMAND}
    -DTEST_PROG=${binaryname}
    -DSOURCEDIR=${CMAKE_CURRENT_SOURCE_DIR}
    -DTEST_SCRIPT=twiss1.madx
    -DTEST_OUTPUT=sample_optics.tfs
    -P testbase.cmake)


add_test(NAME ptc_twiss_1
    COMMAND ${CMAKE_COMMAND}
    -DTEST_PROG=${binaryname}
    -DSOURCEDIR=${CMAKE_CURRENT_SOURCE_DIR}
    -DTEST_SCRIPT=ptc_twiss1.madx
    -DTEST_OUTPUT=twiss_ptc_line.tfs
    -P testbase.cmake)


add_test(NAME ptc_normal_1
    COMMAND ${CMAKE_COMMAND}
    -DTEST_PROG=${binaryname}
    -DSOURCEDIR=${CMAKE_CURRENT_SOURCE_DIR}
    -DTEST_SCRIPT=test_5frs.madx
    -DTEST_OUTPUT=ptc_normal_results.tfs
    -P testbase.cmake)

add_test(NAME ibs_1
    COMMAND ${CMAKE_COMMAND}
    -DTEST_PROG=${binaryname}
    -DSOURCEDIR=${CMAKE_CURRENT_SOURCE_DIR}
    -DTEST_SCRIPT=ibs.mad
    -DTEST_OUTPUT=ibs_output.tfs
    -P testbase.cmake)
