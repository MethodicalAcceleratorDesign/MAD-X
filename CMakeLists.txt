cmake_minimum_required(VERSION 2.6)

PROJECT(madX Fortran C CXX) 
#add name of languages used in the project as well

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
#   cmake_policy(SET CMP0014 OLD) #ignore not having CMakeLists.txt in subfolders (lib/lib64)
endif(COMMAND cmake_policy)


# project version
SET( PROJECT_MAJOR_VERSION 5 )
SET( PROJECT_MINOR_VERSION 00 )
SET( PROJECT_PATCH_LEVEL 05 )

# add source dir as a place for CMake modules (e.g. FindGSL.cmake)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmakesrc") 


# project options
if( NOT BUILD_SHARED_LIBS )
	OPTION( BUILD_SHARED_LIBS "Turn on to build dynamic libraries"  OFF )
endif( NOT BUILD_SHARED_LIBS )

if (APPLE OR BUILD_SHARED_LIBS)
	option( MADX_STATIC "Turn on for static linking" OFF)
else (APPLE OR BUILD_SHARED_LIBS)
	OPTION ( MADX_STATIC "Turn on for static linking" ON)
endif (APPLE OR BUILD_SHARED_LIBS)

#Mad-X specific options (arch. specific options can be added in similar manner):
OPTION ( MADX_NTPSA "Build with NTPSA" ON)
OPTION ( MADX_FORCE_32 "Force 32bit build" OFF )
OPTION ( MADX_FEDORA_FIX "Fix for Fedora>11 for ifort compiler" OFF )

OPTION ( MADX_GOTOBLAS2 "Build with the GOTOBLAS2 libraries" OFF )

if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
	OPTION ( MADX_ONLINE "Build with Online model" ON)
else (CMAKE_SYSTEM_NAME STREQUAL "Linux")
	OPTION ( MADX_ONLINE "Build with Online model" OFF)
endif (CMAKE_SYSTEM_NAME STREQUAL "Linux")

# Default build type (defines different sets of flags)
IF(NOT CMAKE_BUILD_TYPE)
				SET(CMAKE_BUILD_TYPE Release CACHE STRING
								"Choose the type of build, options are: Debug Release"
								FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)



if ( NOT PROJECT_PATCH_LEVEL EQUAL 00 )
 message( "Building a development version" )
 SET (BINARY_POSTFIX "_dev")
 SET (PKG_POSTFIX "-dev")
endif ( NOT PROJECT_PATCH_LEVEL EQUAL 00 )



#adding flags to our compilers:

# 
# we call a subscript that sets up our compiler specific flags...
# 
include(setupCompilerSpecifics)


# General compile flags:
set (CMAKE_C_FLAGS_DEBUG   " ${CMAKE_C_FLAGS_DEBUG} -Wall -pedantic ")
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -funroll-loops -D_CATCH_MEM -D_WRAP_FORTRAN_CALLS -D_WRAP_C_CALLS -D_FULL -I. -I${CMAKE_CURRENT_SOURCE_DIR} ")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -funroll-loops -D_CATCH_MEM -D_WRAP_FORTRAN_CALLS -D_WRAP_C_CALLS -D_FULL -I. -I${CMAKE_CURRENT_SOURCE_DIR} ") #needed for c++ linking
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -D_CATCH_MEM -D_WRAP_FORTRAN_CALLS -D_WRAP_C_CALLS -D_FULL -I. -I${CMAKE_CURRENT_SOURCE_DIR} ") 

# we call a subscript that will set up additional flags depending on 
# architecture chosen.
include(setupArchSpecifics)

if (MADX_ONLINE )
  message("Online Model turned on" )
	if ( NOT MADX_STATIC )
		message( WARNING "You might have problems finding the shared libraries for SDDS" )
	endif ( NOT MADX_STATIC )
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_ONLINE ")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -D_ONLINE ")
endif (MADX_ONLINE )

# OSX specifics:
if (APPLE)
  set(CMAKE_LIBRARY_PATH /usr/lib/ /usr/X11/lib/ ${CMAKE_LIBRARY_PATH})
	set(CMAKE_OSX_DEPLOYMENT_TARGET 10.5)
endif (APPLE)


# C stuff:
if (CMAKE_C_COMPILER_ID MATCHES "GNU" AND NOT CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
# this will probably crash on windows...
  execute_process(COMMAND ${C_COMPILER_NAME} -print-search-dirs
                OUTPUT_VARIABLE gccsearchdirs)
  string(REGEX REPLACE ".*libraries: =(.*)\n"  "\\1" gcclibs "${gccsearchdirs}")
  # need to do this many times because lf95 segfaults on lists with :
  string(REPLACE "/:/"  "/ -L/" gcclibs "${gcclibs}")
  # adding these to the linking process which is handled by a non-gnu fortran compiler in your case
  LINK_DIRECTORIES(${gcclibs}) 
endif (CMAKE_C_COMPILER_ID MATCHES "GNU" AND NOT CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
# end C stuff

#
# we call a subscript that sets up the source files:
# 
include(setupSources)


if(APPLE)
	SET(MACOSX_BUNDLE_STARTUP_COMMAND madx${BINARY_POSTFIX})
  SET(MACOSX_BUNDLE_ICON_FILE "${CMAKE_CURRENT_SOURCE_DIR}/cmakesrc/MadX.icns")
  SET(MACOSX_BUNDLE_LONG_VERSION_STRING "MadX  version ${madX_MAJOR_VERSION}.${madX_MINOR_VERSION}.${madX_PATCH_LEVEL}")
  SET(MACOSX_BUNDLE_BUNDLE_NAME "MadX${BINARY_POSTFIX}")
  SET(MACOSX_BUNDLE_GUI_IDENTIFIER "MadX${BINARY_POSTFIX}")
  # add icns to the .app/Resources with these two commands:
  SET(srcfiles ${srcfiles} ${CMAKE_CURRENT_SOURCE_DIR}/cmakesrc/MadX.icns)
  SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_SOURCE_DIR}/cmakesrc/MadX.icns PROPERTIES MACOSX_PACKAGE_LOCATION Resources)
endif(APPLE)


# adding library:
ADD_LIBRARY(madx ${srcfiles})
# not sure if this is needed...
SET_TARGET_PROPERTIES(madx PROPERTIES LINKER_LANGUAGE Fortran) 

# adding an executable:
add_executable(madxbin MACOSX_BUNDLE  madx_main.F90)
SET_TARGET_PROPERTIES(madxbin PROPERTIES LINKER_LANGUAGE Fortran)

# we want to add _dev to the binary in case this is a dev version for the automatic packaging...
SET_TARGET_PROPERTIES(madxbin PROPERTIES OUTPUT_NAME "madx${BINARY_POSTFIX}")
get_target_property(binaryname madxbin OUTPUT_NAME)

# we need to link executable to our own library:
target_link_libraries(madxbin madx)

# I turn off search for libraries in case you are on Linux,
# to make sure we make use of the lib/lib64 folders
if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux")
#necessary to search for X11 for OSX instead of directly including it
FIND_PACKAGE(X11)
IF(X11_FOUND)
   message("Found X11 libraries")
   INCLUDE_DIRECTORIES(${X11_INCLUDE_DIR})
	 TARGET_LINK_LIBRARIES(madx ${X11_X11_LIB})
ENDIF(X11_FOUND)
else (NOT CMAKE_SYSTEM_NAME STREQUAL "Linux")
  target_link_libraries(madx X11)
endif (NOT CMAKE_SYSTEM_NAME STREQUAL "Linux")

if(LAPACK_FOUND)
    target_link_libraries(madx ${LAPACK_LIBRARIES})
endif()

# Online libraries:
if (MADX_ONLINE )
	target_link_libraries(madx SDDS1c  SDDS1 rpnlib mdbmth mdblib gsl) 
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
	target_link_libraries(madx imf)
endif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
endif (MADX_ONLINE)
target_link_libraries(madx z)
# new additions from Frank:
# On Debian Stable, stdc++ is not linked by default...
target_link_libraries(madx pthread c stdc++ gcc_eh)


# 
# we call a subscript that sets up installation stuff:
# 
include(setupInstallation)


# Testing:
# enable dashboard scripting

# Define useful build name for cdash:
find_program(UNAME NAMES uname)
macro(getuname name flag)
  exec_program("${UNAME}" ARGS "${flag}" OUTPUT_VARIABLE "${name}")
endmacro(getuname)

getuname(osname -s)
# getuname(osrel  -r)
getuname(cpu    -m)

if(CMAKE_Fortran_COMPILER_ID MATCHES "unknown")
    STRING(REGEX REPLACE ".*/([^ ])" "\\1" fcompilername "${CMAKE_Fortran_COMPILER}" )
    SET(BUILDNAME "${osname}-${cpu}-${fcompilername}")
else(CMAKE_Fortran_COMPILER_ID MATCHES "unknown")
    SET(BUILDNAME "${osname}-${cpu}-${CMAKE_Fortran_COMPILER_ID}")
endif(CMAKE_Fortran_COMPILER_ID MATCHES "unknown")

if (MADX_FORCE_32 AND NOT ${CMAKE_SIZEOF_VOID_P} EQUAL 4)
    set(BUILDNAME "${BUILDNAME}-32")
endif (MADX_FORCE_32 AND NOT ${CMAKE_SIZEOF_VOID_P} EQUAL 4)
IF(BUILD_TESTING)
    SET(BUILDNAME "${BUILDNAME}" CACHE STRING "Name of build on the dashboard")
    MARK_AS_ADVANCED(BUILDNAME)
ENDIF(BUILD_TESTING)
# End of defining build name...

enable_testing()
include(CTest)
add_subdirectory(cmakesrc/ctests)
