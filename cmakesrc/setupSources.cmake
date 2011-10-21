
# list of source files
file(GLOB srcfiles *.c *.f90 *.F90 )

file(GLOB main_file mad_main.c)
if(WIN32 OR CYGWIN)
    file(GLOB to_remove gxx11c.c gxx11.f90)
else()
    file(GLOB to_remove gxx11psc.c gxx11ps.f90)
endif()

list(REMOVE_ITEM srcfiles ${main_file} ${to_remove})

message("sources: ${srcfiles}")
if(NOT MADX_FORCE_32)
    # find laplas and blas
    find_package(LAPACK)
endif()

if(LAPACK_FOUND AND BLAS_FOUND)
    message(STATUS "LAPACK uses ${LAPACK_LIBRARIES}")
    file(GLOB to_remove matchlib.f90 matchlib2.f90)
    list(REMOVE_ITEM srcfiles ${to_remove})
else()
    # Note, this only APPENDS -O0 to the compile flags, same as Makefile currently does.
    set_source_files_properties(matchlib2.f90 PROPERTIES COMPILE_FLAGS "-O0")
endif()

# remove source files according to NTPSA option...
if (MADX_NTPSA )
  message(STATUS "NTPSA turned on")
  file(GLOB to_remove c_dabnew.f90)
  set(csrcfiles ${csrcfiles} tpsa.cpp)
else (MADX_NTPSA )
  file(GLOB to_remove c_dabnew_berz.f90 c_tpsa_interface.F90)
endif  (MADX_NTPSA )
list(REMOVE_ITEM srcfiles ${to_remove})

# header files...
file(GLOB headerfiles *.h)
file(GLOB to_remove doxygen.h)
list(REMOVE_ITEM headerfiles ${to_remove})

