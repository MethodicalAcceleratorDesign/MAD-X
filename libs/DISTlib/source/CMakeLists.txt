file(GLOB src_files *.c *.h)

add_library(DISTlib STATIC ${src_files})
if (NOT BUILD_SHARED_LIBS)
    install(TARGETS DISTlib ARCHIVE DESTINATION lib)
endif ()
