
# list of c source files
set(csrcfiles madxp.c matchptcknobs.c rplot.c )
if(WIN32 OR CYGWIN)
    set(csrcfiles ${csrcfiles} gxx11psc.c )
else()
    set(csrcfiles ${csrcfiles} gxx11c.c )
endif()
# list of fortran source files
file(GLOB fsrcfiles *.f90 *.F90)

# Remove some files which should not be compiled..
if(WIN32 OR CYGWIN)
    file(GLOB to_remove gxx11.f90)
else()
    file(GLOB to_remove gxx11ps.f90)
endif()
list(REMOVE_ITEM fsrcfiles ${to_remove})

if(NOT MADX_FORCE_32)
    # find laplas and blas
    find_package(LAPACK)
endif()

if(LAPACK_FOUND AND BLAS_FOUND)
    message(STATUS "LAPACK uses ${LAPACK_LIBRARIES}")
    file(GLOB to_remove matchlib.f90 matchlib2.f90)
    list(REMOVE_ITEM fsrcfiles ${to_remove})
else()
    # Note, this only APPENDS -O0 to the compile flags, same as Makefile currently does.
    set_source_files_properties(matchlib2.f90 PROPERTIES COMPILE_FLAGS "-O0")
endif()

# add source files according to NTPSA option...
if (MADX_NTPSA )
  message(STATUS "NTPSA turned on")
  file(GLOB to_remove c_dabnew.f90)
  set(csrcfiles ${csrcfiles} tpsa.cpp)
else (MADX_NTPSA )
  file(GLOB to_remove c_dabnew_berz.f90 c_tpsa_interface.F90)
endif  (MADX_NTPSA )
list(REMOVE_ITEM fsrcfiles ${to_remove})

include(setupWrappers)

# We need to include the wrapper headers..
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# main source files... 
# TODO: please check this list at some point!
set(srcfiles ${csrcfiles} ${fsrcfiles} ${CMAKE_CURRENT_BINARY_DIR}/fortran_wrappers.c ${CMAKE_CURRENT_BINARY_DIR}/c_wrappers.c )

# header files...
file(GLOB headerfiles *.h)
file(GLOB to_remove doxygen.h)
list(REMOVE_ITEM headerfiles ${to_remove})

