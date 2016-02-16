
# List output files..
file(GLOB TEST_OUTPUT RELATIVE ${SOURCEDIR} "${SOURCEDIR}/*.ref.gz" "${SOURCEDIR}/*.ref")
string(REGEX REPLACE ".ref(.gz)?" ".out" OLD_OUTPUT "${TEST_OUTPUT}")
string(REGEX REPLACE ".ref(.gz)?" "" TEST_OUTPUT "${TEST_OUTPUT}")


# String version for message:
string(REGEX REPLACE ";" " " STRING_TEST_OUTPUT "${TEST_OUTPUT}")

file(REMOVE ${OLD_OUTPUT} ${TEST_OUTPUT})

# Run simulation..
message("COMMAND ${TEST_PROG} < ${SOURCEDIR}/${TEST_NAME}.madx")
execute_process(COMMAND ${TEST_PROG} INPUT_FILE ${SOURCEDIR}/${TEST_NAME}.madx
   OUTPUT_FILE ${TEST_NAME}.out WORKING_DIRECTORY ${SOURCEDIR} RESULT_VARIABLE HAD_ERROR)

if(HAD_ERROR)
    file(READ ${SOURCEDIR}/${TEST_NAME}.out madxoutput)
    message("${madxoutput}")
    message(FATAL_ERROR "\n\t -- CTEST ERROR -- Test failed with error ${HAD_ERROR}")
else()
    # Run ndiff if all went well..
    message("COMMAND ${NDIFF} -b -l -t ${TEST_NAME} ${STRING_TEST_OUTPUT}")
    execute_process(COMMAND ${NDIFF} -b -l -t ${TEST_NAME} ${TEST_OUTPUT}
        WORKING_DIRECTORY ${SOURCEDIR} RESULT_VARIABLE NDIFF_ERROR)
    if(NDIFF_ERROR)
       message(FATAL_ERROR "Test failed with ndiff error ${NDIFF_ERROR}")
    endif()
endif()

