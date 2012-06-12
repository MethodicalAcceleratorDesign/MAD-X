execute_process(COMMAND cat ${SOURCEDIR}/${TEST_NAME}.madx COMMAND ${TEST_PROG}
   OUTPUT_FILE "out.txt" WORKING_DIRECTORY ${SOURCEDIR} RESULT_VARIABLE HAD_ERROR)

if(HAD_ERROR)
    message(FATAL_ERROR "Test failed with error ${HAD_ERROR}")
 else()
    execute_process(COMMAND ${NUMDIFF} -q -b -c -n -t ${TEST_NAME} out.txt ref.txt cfg.txt  WORKING_DIRECTORY ${SOURCEDIR} RESULT_VARIABLE NUMDIFF_ERROR)
    if(NUMDIFF_ERROR)
       message(FATAL_ERROR "Test failed with numdiff error ${NUMDIFF_ERROR}")
    endif()
endif()


