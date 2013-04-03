
# List output files..
file(GLOB TEST_OUTPUT RELATIVE ${SOURCEDIR} ${SOURCEDIR}/*.ref)
string(REGEX REPLACE ".ref" ".out" OLD_OUTPUT "${TEST_OUTPUT}")
string(REGEX REPLACE ".ref" "" TEST_OUTPUT "${TEST_OUTPUT}")

# String version for message:
string(REGEX REPLACE ";" " " STRING_TEST_OUTPUT "${TEST_OUTPUT}")

file(REMOVE ${OLD_OUTPUT} ${TEST_OUTPUT})

# Run simulation..
message("COMMAND ${TEST_PROG} < ${SOURCEDIR}/${TEST_NAME}.madx")
execute_process(COMMAND ${TEST_PROG} INPUT_FILE ${SOURCEDIR}/${TEST_NAME}.madx
   OUTPUT_FILE ${TEST_NAME}.out WORKING_DIRECTORY ${SOURCEDIR} RESULT_VARIABLE HAD_ERROR)

if(HAD_ERROR)
    message(FATAL_ERROR "Test failed with error ${HAD_ERROR}")
else()
    # Run numdiff if all went well..
    message("COMMAND ${NUMDIFF} -q -b --xcheck -l -n -t ${TEST_NAME} ${STRING_TEST_OUTPUT}")
    execute_process(COMMAND ${NUMDIFF} -q -b --xcheck -l -n -t ${TEST_NAME} ${TEST_OUTPUT}
        WORKING_DIRECTORY ${SOURCEDIR} RESULT_VARIABLE NUMDIFF_ERROR)
    if(NUMDIFF_ERROR)
       message(FATAL_ERROR "Test failed with numdiff error ${NUMDIFF_ERROR}")
    endif()
endif()

