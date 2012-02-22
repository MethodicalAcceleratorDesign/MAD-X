
# list of source files
file(GLOB mad_core_src mad_*.c *.f90 *.F90 *.cpp )
file(GLOB mad_src 
          S[p-t]_*.f90 
          Sra_*.f90 
          madx_ptc*.f90 
          ptc_*.f90 
          user2_photon.f90 
          wrap.f90 
          mad_ptc*.c 
          mad_[cmest]*.c)
file(GLOB mad_fpp_src [a-z]_*.f90 [a-z]_*.F90 tpsa.cpp)
file(GLOB mad_ptc_core_src S[a-f]_*.f90)
file(GLOB mad_ptc_src S[g-o]_*.f90 Sma*.f90)

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
else (MADX_NTPSA )
   file(GLOB ntpsa_remove c_dabnew_berz.f90 c_tpsa_interface.f90 tpsa.cpp)
endif  (MADX_NTPSA )

set(remove_list ${main_file} ${gxx_remove} ${sdds_remove} ${ntpsa_remove} ${lapack_remove})
list(REMOVE_ITEM mad_core_src ${remove_list} ${mad_src} ${mad_fpp_src} ${mad_ptc_core_src} ${mad_ptc_src})
list(REMOVE_ITEM mad_src ${remove_list})
list(REMOVE_ITEM mad_fpp_src ${remove_list})
list(REMOVE_ITEM mad_ptc_core_src ${remove_list})
list(REMOVE_ITEM mad_ptc_src ${remove_list})

# header files...
file(GLOB headerfiles *.h)

# adding library:
add_library(madx ${mad_src})
add_library(madx_core ${mad_core_src})
add_library(madx_fpp ${mad_fpp_src})
add_library(madx_ptc_core ${mad_ptc_core_src})
add_library(madx_ptc ${mad_ptc_src})

target_link_libraries(madx_ptc_core madx_fpp)
target_link_libraries(madx_ptc madx_ptc_core)
target_link_libraries(madx madx_ptc)
target_link_libraries(madx madx_core)


# not sure if this is needed...
set_target_properties(madx PROPERTIES LINKER_LANGUAGE Fortran)

# set library version:
set_target_properties(madx     PROPERTIES VERSION ${MADX_VERSION} SOVERSION 1)
set_target_properties(madx_fpp PROPERTIES VERSION ${MADX_VERSION} SOVERSION 1)

# adding an executable:
if(MADX_BUNDLE)
   add_executable(madxbin MACOSX_BUNDLE  mad_main.c)
else()
   add_executable(madxbin mad_main.c)
endif()
set_target_properties(madxbin PROPERTIES LINKER_LANGUAGE Fortran VERSION ${MADX_VERSION})

# we want to add _dev to the binary in case this is a dev version for the automatic packaging...
set_target_properties(madxbin PROPERTIES OUTPUT_NAME "madx${BINARY_POSTFIX}")

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
