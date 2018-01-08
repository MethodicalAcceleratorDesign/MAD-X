

file(GLOB src_files *.f90 *.F90 *.cpp)

# remove source files according to NTPSA option...
if(MADX_NTPSA)
   message(STATUS "NTPSA turned on")
   file(GLOB ntpsa_remove c_dabnew.f90)
else()
   file(GLOB ntpsa_remove c_dabnew_berz.f90 c_tpsa_interface.F90 tpsa.cpp)
endif()
list(REMOVE_ITEM src_files ${ntpsa_remove})


# Some weird problem with this file with gfortran < 4.6..
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND 
        CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "4.6")
    set_source_files_properties(l_complex_taylor.f90 COMPILE_FLAGS "-O2")
endif()


add_library(ptc ${src_files})

set_target_properties(ptc PROPERTIES VERSION ${MADX_VERSION} SOVERSION 1)

madx_install_targets(ptc)
