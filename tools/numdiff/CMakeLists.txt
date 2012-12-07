cmake_minimum_required(VERSION 2.4)

PROJECT(numdiff C)

if(NOT "${CMAKE_C_FLAGS}" MATCHES "-std=c99")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
endif()

file(GLOB numdiff_src src/*.c)

include_directories(src)
add_executable(numdiff ${numdiff_src})
target_link_libraries(numdiff m)


install(TARGETS numdiff
  RUNTIME
  DESTINATION bin
  COMPONENT Runtime
  )
