
# We want to run coverage analysis of
# the files compiled here..
include(buildtypeCoverage)

# list of source files
file(GLOB src_files
   mad_*.c
   *.f90
   *.F90
   *.cpp)
file(GLOB main_file mad_main.c)
if(WIN32 OR CYGWIN)
   file(GLOB gxx_remove gxx11c.c gxx11.f90)
else()
   file(GLOB gxx_remove gxx11psc.c gxx11ps.f90)
endif()

# Some source files which should not be included:
file(GLOB sdds_remove matchsdds.c)
# the old makethin should not be used:
file(GLOB mkthin_remove mad_mkthin.c)

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

set(remove_list ${gxx_remove} ${sdds_remove} ${mkthin_remove} ${lapack_remove})
list(REMOVE_ITEM src_files ${remove_list})

# header files...
file(GLOB headerfiles *.h)

# adding library:
add_library(madx ${src_files})

target_link_libraries(madx ptc)

# Garbage Collector:
if(USE_GC)
   target_link_libraries(madx gc-lib)
endif()

# not sure if this is needed...
set_target_properties(madx PROPERTIES LINKER_LANGUAGE Fortran)

# set library version:
set_target_properties(madx     PROPERTIES VERSION ${MADX_VERSION})
set_target_properties(madx     PROPERTIES CXX_STANDARD 11)


# adding an executable:
if(MADX_BUNDLE)
   add_executable(madxbin MACOSX_BUNDLE  ${main_file})
else()
   add_executable(madxbin ${main_file})
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
if(NOT WIN32)
   target_link_libraries(madx z)
endif()

find_package(Threads)
find_package(STDCPP)
target_link_libraries(madx  ${CMAKE_THREAD_LIBS_INIT} ${STDCPP_LIBS})

# Garbage collector:
if(USE_GC)
   include_directories(${CMAKE_SOURCE_DIR}/libs/gc/gc-7.2f/include)
endif()

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
