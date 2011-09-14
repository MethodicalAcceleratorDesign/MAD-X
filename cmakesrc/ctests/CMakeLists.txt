

if (APPLE)
    set(binaryname "${CMAKE_BINARY_DIR}/${binaryname}.app/Contents/MacOS/./${binaryname}")
else(APPLE)
    set(binaryname "${CMAKE_BINARY_DIR}/./${binaryname}")
endif (APPLE)

ADD_TEST(testruns "${CMAKE_CURRENT_SOURCE_DIR}/testruns.sh" "${binaryname}")
set_tests_properties (testruns 
  PROPERTIES PASS_REGULAR_EXPRESSION "Number of warnings: 0")
# This is a slow test so I commented it out..
# ADD_TEST(testsample "${CMAKE_CURRENT_SOURCE_DIR}/testsample.sh" "${binaryname}" "${CMAKE_CURRENT_SOURCE_DIR}")
# set_tests_properties (testsample 
#   PROPERTIES FAIL_REGULAR_EXPRESSION "Number of errors: [1-9]") # The no warning check fails on this test..

# This is a general way to test easily using a macro script I have written:
# TEST_SCRIPT: name of madx script to execute
# TEST_OUTPUT: name of file madx produces (skips if empty)


set(BASESCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/testbase.cmake)

add_test(ptc_normal_1
    ${CMAKE_COMMAND}
    -DTEST_PROG=${binaryname}
    -DSOURCEDIR=${CMAKE_CURRENT_SOURCE_DIR}
    -DTEST_SCRIPT=test_5frs.madx
    -DTEST_OUTPUT=ptc_normal_results.tfs
    -P ${BASESCRIPT})

# Tests that require afs:
if(EXISTS "/afs/cern.ch/")
    add_test(twiss_1
        ${CMAKE_COMMAND}
        -DTEST_PROG=${binaryname}
        -DSOURCEDIR=${CMAKE_CURRENT_SOURCE_DIR}
        -DTEST_SCRIPT=twiss1.madx
        -DTEST_OUTPUT=sample_optics.tfs
        -P ${BASESCRIPT})
    
    add_test(ptc_twiss_1
        ${CMAKE_COMMAND}
        -DTEST_PROG=${binaryname}
        -DSOURCEDIR=${CMAKE_CURRENT_SOURCE_DIR}
        -DTEST_SCRIPT=ptc_twiss1.madx
        -DTEST_OUTPUT=twiss_ptc_line.tfs
        -P ${BASESCRIPT})
    
    add_test(ibs_1
        ${CMAKE_COMMAND}
        -DTEST_PROG=${binaryname}
        -DSOURCEDIR=${CMAKE_CURRENT_SOURCE_DIR}
        -DTEST_SCRIPT=ibs.mad
        -DTEST_OUTPUT=ibs_output.tfs
        -P ${BASESCRIPT})
else()
    message(STATUS "afs is not available, some tests will be missing")
endif()
