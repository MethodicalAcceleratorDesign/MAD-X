
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

# adding library:
add_library(madx ${srcfiles})
# not sure if this is needed...
set_target_properties(madx PROPERTIES LINKER_LANGUAGE Fortran) 

# adding an executable:
add_executable(madxbin MACOSX_BUNDLE  mad_main.c)
set_target_properties(madxbin PROPERTIES LINKER_LANGUAGE Fortran)

# we want to add _dev to the binary in case this is a dev version for the automatic packaging...
set_target_properties(madxbin PROPERTIES OUTPUT_NAME "madx${BINARY_POSTFIX}")
get_target_property(binaryname madxbin OUTPUT_NAME)

# we need to link executable to our own library:
target_link_libraries(madxbin madx)

# I turn off search for libraries in case you are on Linux,
# to make sure we make use of the lib/lib64 folders
if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux")
    #necessary to search for X11 for OSX instead of directly including it
    find_package(X11)
    if(X11_FOUND)
        include_directories(${X11_INCLUDE_DIR})
        target_link_libraries(madx ${X11_X11_LIB})
    endif()
else()
    target_link_libraries(madx X11)
endif()

if(LAPACK_FOUND AND BLAS_FOUND)
    target_link_libraries(madx ${LAPACK_LIBRARIES})
endif()

# Online libraries:
if(MADX_ONLINE)
    target_link_libraries(madx ${SDDS_LIBRARIES}) 
endif()
target_link_libraries(madx z)
# new additions from Frank:
# On Debian Stable, stdc++ is not linked by default...
target_link_libraries(madx pthread c stdc++ gcc_eh)

# Defining specific stuff for packaging:
# Should be called after sources are defined,
# before targets are defined.
include(setupPackaging)

# 
# we call a subscript that sets up installation stuff:
# 
# Should be called after targets are defined
# 
include(setupInstallation)
