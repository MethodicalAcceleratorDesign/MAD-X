execute_process(COMMAND cat ${SOURCEDIR}/${TEST_SCRIPT} COMMAND ${TEST_PROG}
                RESULT_VARIABLE HAD_ERROR)
if(HAD_ERROR)
    message(FATAL_ERROR "Test failed with error ${HAD_ERROR}")
endif()

# This will always fail because date/version is included in output...
if(TEST_OUTPUT)
    find_package(PythonInterp)
    if(PYTHONINTERP_FOUND)
        execute_process(COMMAND ${PYTHON_EXECUTABLE} ${SOURCEDIR}/checkConsistency.py
            ${SOURCEDIR}/madx_testing_output/${TEST_OUTPUT} ${TEST_OUTPUT}
            RESULT_VARIABLE DIFFERENT)
        if(DIFFERENT)
            message(FATAL_ERROR "Test failed - files differ. Got return value ${DIFFERENT}")
        endif()
    endif()
endif()
