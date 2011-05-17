cmake_minimum_required(VERSION 2.6)

PROJECT(madX C CXX Fortran) 
#add name of languages used in the project as well

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
#   cmake_policy(SET CMP0014 OLD) #ignore not having CMakeLists.txt in subfolders (lib/lib64)
endif(COMMAND cmake_policy)


# project version
SET( PROJECT_MAJOR_VERSION 5 )
SET( PROJECT_MINOR_VERSION 00 )
SET( PROJECT_PATCH_LEVEL 00 )

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
get_filename_component (C_COMPILER_NAME ${CMAKE_C_COMPILER} NAME)
if (C_COMPILER_NAME MATCHES "gcc" AND NOT CMAKE_Fortran_COMPILER MATCHES "gfortran")
# this will probably crash on windows...
  execute_process(COMMAND ${C_COMPILER_NAME} -print-search-dirs
                OUTPUT_VARIABLE gccsearchdirs)
  string(REGEX REPLACE ".*libraries: =(.*)\n"  "\\1" gcclibs "${gccsearchdirs}")
  # need to do this many times because lf95 segfaults on lists with :
  string(REPLACE "/:/"  "/ -L/" gcclibs "${gcclibs}")
  # adding these to the linking process which is handled by a non-gnu fortran compiler in your case
  LINK_DIRECTORIES(${gcclibs}) 
endif (C_COMPILER_NAME MATCHES "gcc" AND NOT CMAKE_Fortran_COMPILER MATCHES "gfortran")
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
add_executable(madxbin MACOSX_BUNDLE  madxp.c madx_main.F90)
SET_TARGET_PROPERTIES(madxbin PROPERTIES LINKER_LANGUAGE Fortran)

# we want to add _dev to the binary in case this is a dev version for the automatic packaging...
SET_TARGET_PROPERTIES(madxbin PROPERTIES OUTPUT_NAME "madx${BINARY_POSTFIX}")

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

enable_testing()
include(CTest)
set (CTEST_PROJECT_NAME "MadX")
ADD_TEST(testruns "${CMAKE_CURRENT_SOURCE_DIR}/cmakesrc/ctests/testruns.sh" "${BINARY_POSTFIX}")
set_tests_properties (testruns 
  PROPERTIES PASS_REGULAR_EXPRESSION "Number of warnings: 0")
ADD_TEST(testsample "${CMAKE_CURRENT_SOURCE_DIR}/cmakesrc/ctests/testsample.sh" "${BINARY_POSTFIX}" "${CMAKE_CURRENT_SOURCE_DIR}/cmakesrc/ctests")
set_tests_properties (testsample 
  PROPERTIES FAIL_REGULAR_EXPRESSION "Number of errors: [1-9]") # The no warning check fails on this test..

