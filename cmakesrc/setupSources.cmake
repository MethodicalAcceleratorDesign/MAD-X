
# list of c source files
set(csrcfiles madxp.c gxx11c.c matchptcknobs.c rplot.c )
# list of fortran source files
file(GLOB fsrcfiles *.f90 *.F90)

# Remove some files which should never be compiled..
file(GLOB to_remove gxx11ps.f90)
list(REMOVE_ITEM fsrcfiles ${to_remove})

if(MADX_STATIC)
    if(WIN32)
        set(CMAKE_FIND_LIBRARY_SUFFIXES .lib)
    elseif(APPLE)
        set(CMAKE_FIND_LIBRARY_SUFFIXES .a .lib)
    else()
        set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
    endif()
endif()
find_package(LAPACK) # (lapack requires blas...)
if(MADX_RICCARDO_FIX AND NOT LAPACK_FOUND)
    find_library(BLAS_LIBRARIES blas)
    find_library(LAPACK_LIBRARIES lapack)
    if(LAPACK_LIBRARIES AND BLAS_LIBRARIES)
        set(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
        set(LAPACK_FOUND TRUE)
    endif()
endif()
if(LAPACK_FOUND AND BLAS_FOUND)
    file(GLOB to_remove matchlib.f90 matchlib2.f90)
    list(REMOVE_ITEM fsrcfiles ${to_remove})
else()
    # Note, this only APPENDS -O0 to the compile flags, same as Makefile currently does.
    set_source_files_properties(matchlib2.f90 PROPERTIES COMPILE_FLAGS "-O0") 
endif()

# add source files according to NTPSA option...
if (MADX_NTPSA )
  message("NTPSA turned on")
  file(GLOB to_remove c_dabnew.f90)
  set(csrcfiles ${csrcfiles} tpsa.cpp)
else (MADX_NTPSA )
  file(GLOB to_remove c_dabnew_berz.f90 c_tpsa_interface.F90)
endif  (MADX_NTPSA )
list(REMOVE_ITEM fsrcfiles ${to_remove})

#execute python wrapper scripts (you need to be dependent on one of the output files or else this command will never be ran):
# Unsure about dependencies.. Might be an overkill this one.
ADD_CUSTOM_COMMAND(
  OUTPUT c_wrappers.c c_wrappers.h c_prototypes.h c_wrappers_prototypes.h
  DEPENDS ${csrcfiles} ${fsrcfiles}
  COMMAND python2 wrap_C_calls.py -o ${CMAKE_CURRENT_BINARY_DIR} 
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "Creating C wrapper files"
  )
ADD_CUSTOM_COMMAND(
  OUTPUT fortran_wrappers.c fortran_wrappers.h fortran_prototypes.h fortran_wrappers_prototypes.h
  DEPENDS ${csrcfiles} ${fsrcfiles}
  COMMAND python2 wrap_fortran_calls.py --outdir=${CMAKE_CURRENT_BINARY_DIR} 
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "Creating fortran wrapper files"
  )
# execute python wrap_fortran_calls.py
# COMMAND python wrap_C_calls.py


# main source files... 
# TODO: please check this list at some point!
set(srcfiles ${csrcfiles} ${fsrcfiles} ${CMAKE_CURRENT_BINARY_DIR}/fortran_wrappers.c ${CMAKE_CURRENT_BINARY_DIR}/c_wrappers.c )

# header files...
file(GLOB headerfiles *.h)
file(GLOB to_remove doxygen.h)
list(REMOVE_ITEM headerfiles ${to_remove})

