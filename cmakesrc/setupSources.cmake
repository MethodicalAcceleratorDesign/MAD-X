
# list of source files
file(GLOB srcfiles mad_*.c gxx*.c madxp.c matchptcknobs.c rplot.c *.f90 *.F90 )

file(GLOB main_file mad_main.c)
if(WIN32 OR CYGWIN)
    file(GLOB gxx_remove gxx11c.c gxx11.f90)
else()
    file(GLOB gxx_remove gxx11psc.c gxx11ps.f90)
endif()

# Some source files which should not be included:
file(GLOB sdds_remove matchsdds.c)


if(NOT MADX_FORCE_32)
    # find laplas and blas
    find_package(LAPACK)
endif()

if(LAPACK_FOUND AND BLAS_FOUND)
    message(STATUS "LAPACK uses ${LAPACK_LIBRARIES}")
    file(GLOB lapack_remove matchlib.f90 matchlib2.f90)
else()
    # Note, this only APPENDS -O0 to the compile flags, same as Makefile currently does.
    set_source_files_properties(matchlib2.f90 PROPERTIES COMPILE_FLAGS "-O0")
endif()

# remove source files according to NTPSA option...
if (MADX_NTPSA )
  message(STATUS "NTPSA turned on")
  file(GLOB ntpsa_remove c_dabnew.f90)
  set(srcfiles ${srcfiles} tpsa.cpp)
else (MADX_NTPSA )
  file(GLOB ntpsa_remove c_dabnew_berz.f90 c_tpsa_interface.f90)
endif  (MADX_NTPSA )

list(REMOVE_ITEM srcfiles ${main_file} ${gxx_remove} ${sdds_remove} ${ntpsa_remove} ${lapack_remove})

# header files...
file(GLOB headerfiles *.h)
file(GLOB to_remove doxygen.h)
list(REMOVE_ITEM headerfiles ${to_remove})

