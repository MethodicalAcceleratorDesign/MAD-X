## This file should be placed in the root directory of your project.
## Then modify the CMakeLists.txt file in the root directory of your
## project to incorporate the testing dashboard.
## # The following are required to uses Dart and the Cdash dashboard
##   ENABLE_TESTING()
##   INCLUDE(CTest)
set(CTEST_PROJECT_NAME "Mad-X")
set(CTEST_NIGHTLY_START_TIME "00:00:00 EST")
set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "137.138.26.237")
set(CTEST_DROP_LOCATION "/cdash/submit.php?project=Mad-X")
set(CTEST_DROP_SITE_CDASH TRUE)
set(CTEST_UPDATE_COMMAND "svn")



SET(CTEST_MEMORYCHECK_COMMAND "valgrind")
# SET(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE "${CTEST_SCRIPT_DIRECTORY}/CMakeValgrindSuppressions.supp")
SET(CTEST_MEMORYCHECK_COMMAND_OPTIONS "--gen-suppressions=all --child-silent-after-fork=yes --trace-children=yes --trace-children-skip=${valgrind_skip} --track-origins=yes -q --leak-check=yes --show-reachable=yes --workaround-gcc296-bugs=yes --num-callers=100 -v")


if (APPLE)
    set(binaryname "${binaryname}.app/Contents/MacOS/./${binaryname}")
else(APPLE)
    set(binaryname "./${binaryname}")
endif (APPLE)

ADD_TEST(testruns "${CMAKE_CURRENT_SOURCE_DIR}/cmakesrc/ctests/testruns.sh" "${binaryname}")
set_tests_properties (testruns 
  PROPERTIES PASS_REGULAR_EXPRESSION "Number of warnings: 0")

# This is a slow test so I commented it out..
# ADD_TEST(testsample "${CMAKE_CURRENT_SOURCE_DIR}/cmakesrc/ctests/testsample.sh" "${binaryname}" "${CMAKE_CURRENT_SOURCE_DIR}/cmakesrc/ctests")
# set_tests_properties (testsample 
#   PROPERTIES FAIL_REGULAR_EXPRESSION "Number of errors: [1-9]") # The no warning check fails on this test..

# This is a general way to test easily using a macro script I have written:
# TEST_SCRIPT: name of madx script to execute
# TEST_OUTPUT: name of file madx produces (skips if empty)
add_test(NAME twiss_1
    COMMAND ${CMAKE_COMMAND}
    -DTEST_PROG=${binaryname}
    -DSOURCEDIR=${CMAKE_CURRENT_SOURCE_DIR}/cmakesrc/ctests
    -DTEST_SCRIPT=twiss1.madx
    -DTEST_OUTPUT=sample_optics.tfs
    -P ${CMAKE_CURRENT_SOURCE_DIR}/cmakesrc/ctests/testbase.cmake)


add_test(NAME ptc_twiss_1
    COMMAND ${CMAKE_COMMAND}
    -DTEST_PROG=${binaryname}
    -DSOURCEDIR=${CMAKE_CURRENT_SOURCE_DIR}/cmakesrc/ctests
    -DTEST_SCRIPT=ptc_twiss1.madx
    -DTEST_OUTPUT=twiss_ptc_line.tfs
    -P ${CMAKE_CURRENT_SOURCE_DIR}/cmakesrc/ctests/testbase.cmake)


add_test(NAME ptc_normal_1
    COMMAND ${CMAKE_COMMAND}
    -DTEST_PROG=${binaryname}
    -DSOURCEDIR=${CMAKE_CURRENT_SOURCE_DIR}/cmakesrc/ctests
    -DTEST_SCRIPT=test_5frs.madx
    -DTEST_OUTPUT=ptc_normal_results.tfs
    -P ${CMAKE_CURRENT_SOURCE_DIR}/cmakesrc/ctests/testbase.cmake)

add_test(NAME ibs_1
    COMMAND ${CMAKE_COMMAND}
    -DTEST_PROG=${binaryname}
    -DSOURCEDIR=${CMAKE_CURRENT_SOURCE_DIR}/cmakesrc/ctests
    -DTEST_SCRIPT=ibs.mad
    -DTEST_OUTPUT=ibs_output.tfs
    -P ${CMAKE_CURRENT_SOURCE_DIR}/cmakesrc/ctests/testbase.cmake)