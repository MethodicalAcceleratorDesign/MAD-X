

set(numdiff_src src/constraint.c  src/context.c  src/error.c  src/main.c  src/ndiff.c  src/utest.c  src/utils.c)

include_directories(src)

add_executable(numdiff ${numdiff_src})
target_link_libraries(numdiff m)
