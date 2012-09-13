

file(GLOB src_files *.f90)

# Some weird problem with this file with gfortran < 4.6..
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND 
        CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "4.6")
    set_source_files_properties(l_complex_taylor.f90 COMPILE_FLAGS "-O2")
endif()


add_library(ptc ${src_files})

set_target_properties(ptc PROPERTIES VERSION ${MADX_VERSION} SOVERSION 1)

madx_install_targets(ptc)
