cmake_minimum_required(VERSION 2.4)

PROJECT(numdiff C)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")

file(GLOB numdiff_src src/*.c)

include_directories(src)
add_executable(numdiff ${numdiff_src})
target_link_libraries(numdiff m)


install(TARGETS numdiff
  RUNTIME
  DESTINATION bin
  COMPONENT Runtime
  )
