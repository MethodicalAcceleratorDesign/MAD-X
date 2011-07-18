execute_process(COMMAND cat ${SOURCEDIR}/${TEST_SCRIPT} COMMAND ${TEST_PROG}
                RESULT_VARIABLE HAD_ERROR)
if(HAD_ERROR)
    message(FATAL_ERROR "Test failed")
endif()

# This will always fail because date/version is included in output...
# if(TEST_OUTPUT)
#     execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files
#         ${TEST_OUTPUT} /afs/cern.ch/user/y/ylevinse/scratch1/public/madx_testing_output/${TEST_OUTPUT}
#         RESULT_VARIABLE DIFFERENT)
#     if(DIFFERENT)
#         message(FATAL_ERROR "Test failed - files differ")
#     endif()
# endif()
