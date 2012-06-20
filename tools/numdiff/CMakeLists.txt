
file(GLOB numdiff_src src/*.c)

include_directories(src)

add_executable(numdiff ${numdiff_src})
target_link_libraries(numdiff m)
