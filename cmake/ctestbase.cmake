
# We must make a list so it is understood as multiple
# input arguments..
string(REGEX REPLACE " " ";" _TEST_OUTPUT ${TEST_OUTPUT})

# Remove files in case this test was ran before..
foreach(_f ${_TEST_OUTPUT} ${TEST_NAME}.out)
   file(REMOVE ${SOURCEDIR}/${_f})
endforeach()

# Run simulation..
message("COMMAND ${TEST_PROG} < ${SOURCEDIR}/${TEST_NAME}.madx")
execute_process(COMMAND ${TEST_PROG} INPUT_FILE ${SOURCEDIR}/${TEST_NAME}.madx
   OUTPUT_FILE ${TEST_NAME}.out WORKING_DIRECTORY ${SOURCEDIR} RESULT_VARIABLE HAD_ERROR)

if(HAD_ERROR)
    message(FATAL_ERROR "Test failed with error ${HAD_ERROR}")
else()
    # Run numdiff if all went well..
    message("COMMAND ${NUMDIFF} -q -b -c -l -n -t ${TEST_NAME} ${TEST_OUTPUT}")
    execute_process(COMMAND ${NUMDIFF} -q -b -c -l -n -t ${TEST_NAME} ${_TEST_OUTPUT}
        WORKING_DIRECTORY ${SOURCEDIR} RESULT_VARIABLE NUMDIFF_ERROR)
    if(NUMDIFF_ERROR)
       message(FATAL_ERROR "Test failed with numdiff error ${NUMDIFF_ERROR}")
    endif()
endif()


