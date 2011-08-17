# small script that detects number of processors.
# after this script, you can add the following
# to your ctest script:
#
# if(PROCESSOR_COUNT)
#  set(CTEST_BUILD_FLAGS "-j${PROCESSOR_COUNT}")
# endif()
#
# Source: http://www.kitware.com/blog/home/post/63
if(NOT DEFINED PROCESSOR_COUNT)
  # Unknown:
  set(PROCESSOR_COUNT 0)

  # Linux:
  set(cpuinfo_file "/proc/cpuinfo")
  if(EXISTS "${cpuinfo_file}")
    file(STRINGS "${cpuinfo_file}" procs REGEX "^processor.: [0-9]+$")
    list(LENGTH procs PROCESSOR_COUNT)
  endif()

  # Mac:
  if(APPLE)
    find_program(cmd_sys_pro "system_profiler")
    if(cmd_sys_pro)
      execute_process(COMMAND ${cmd_sys_pro} OUTPUT_VARIABLE info)
      string(REGEX REPLACE "^.*Total Number Of Cores: ([0-9]+).*$" "\\1"
        PROCESSOR_COUNT "${info}")
    endif()
  endif()

  # Windows:
  if(WIN32)
    set(PROCESSOR_COUNT "$ENV{NUMBER_OF_PROCESSORS}")
  endif()
endif()

message(STATUS "Number of processors found: ${PROCESSOR_COUNT}")
