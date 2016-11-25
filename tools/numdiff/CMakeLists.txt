cmake_minimum_required(VERSION 2.4)

project(ndiff C)

option( INSTALL_NDIFF "Add ndiff install target" ON)

if(NOT "${CMAKE_C_FLAGS}" MATCHES "std=c99")
   if(NOT (WIN32 AND CMAKE_C_COMPILER_ID STREQUAL "Intel"))
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
   else()
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /Qstd=c99")
   endif()
endif()

file(GLOB ndiff_src src/*.c)

include_directories(src)
add_executable(ndiff ${ndiff_src})
if(UNIX)
   target_link_libraries(ndiff m)
endif()

if(INSTALL_NDIFF)
   install(TARGETS ndiff
      RUNTIME
      DESTINATION bin
      COMPONENT Runtime
      )
endif()

# Define tests (without using config from parent folder)
enable_testing()
foreach(testname test-register test-register-2 test-register-3)
    set(testpath  "${CMAKE_SOURCE_DIR}/tests/${testname}/${testname}")
    add_test(${testname} ndiff ${testpath}.nd.lhs ${testpath}.nd.rhs ${testpath}.nd.cfg)
    set_tests_properties(${testname} PROPERTIES PASS_REGULAR_EXPRESSION "diffs detected over")
endforeach()
