
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

# find laplas and blas
find_package(LAPACK)

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

#execute python wrapper scripts (you need to be dependent on one of the output files or else this command will never be ran):
# Unsure about dependencies.. Might be an overkill this one.
find_package(PythonInterp REQUIRED)
ADD_CUSTOM_COMMAND(
  OUTPUT c_wrappers.c c_wrappers.h c_prototypes.h c_wrappers_prototypes.h
  DEPENDS ${csrcfiles} ${fsrcfiles}
  COMMAND ${PYTHON_EXECUTABLE} wrap_C_calls.py -o ${CMAKE_CURRENT_BINARY_DIR} 
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "Creating C wrapper files"
  )
ADD_CUSTOM_COMMAND(
  OUTPUT fortran_wrappers.c fortran_wrappers.h fortran_prototypes.h fortran_wrappers_prototypes.h
  DEPENDS ${csrcfiles} ${fsrcfiles}
  COMMAND ${PYTHON_EXECUTABLE} wrap_fortran_calls.py --outdir=${CMAKE_CURRENT_BINARY_DIR} 
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "Creating fortran wrapper files"
  )
# execute python wrap_fortran_calls.py
# COMMAND python wrap_C_calls.py

# We need to include the wrapper headers..
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# main source files... 
# TODO: please check this list at some point!
set(srcfiles ${csrcfiles} ${fsrcfiles} ${CMAKE_CURRENT_BINARY_DIR}/fortran_wrappers.c ${CMAKE_CURRENT_BINARY_DIR}/c_wrappers.c )

# header files...
file(GLOB headerfiles *.h)
file(GLOB to_remove doxygen.h)
list(REMOVE_ITEM headerfiles ${to_remove})

