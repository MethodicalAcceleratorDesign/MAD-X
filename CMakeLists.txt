# Suppress warnings on cygwin
# win32 is not anymore equal to 1 on cygwin
# with newer cmake. We still assume so though..
# Note:
# - In cygwin, the boolean CYGWIN is defined and 1
# - On Windows, WIN32 bool is defined, but 
# after cmake 2.8.4 it's normally zero when using
# cygwin.
set(CMAKE_LEGACY_CYGWIN_WIN32 0)

cmake_minimum_required(VERSION 2.8)

PROJECT(madX C CXX Fortran)
#add name of languages used in the project as well

# Some policies regarding configure
# warnings we do not want..
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
  cmake_policy(SET CMP0042 NEW)
endif(COMMAND cmake_policy)

# add source dir as a place for CMake modules (e.g. FindGSL.cmake)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake") 

# When we don't need to support cmake <2.8.3 anymore, remove these three lines:
if(CMAKE_VERSION VERSION_LESS "2.8.3")
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/compilers") 
endif()

# Some initializations specific for Mad-X
include(madxInit)

# predefined macros...
include(install_macros)

# Defining project options:
include(setupOptions)

# We add libraries before setting up compilers because
# some of our options are incompatible with e.g. GC
if(USE_GC)
   add_subdirectory(libs/gc/gc)
endif()

# 
# we call a subscript that sets up our compiler specific flags...
# 
include(setupCompilerSpecifics)


# we call a subscript that will set up additional flags depending on 
# architecture chosen.
include(setupArchSpecifics)


# Add sources subdirectory:
add_subdirectory(libs/ptc)
add_subdirectory(src)
add_subdirectory(tools)
add_subdirectory(syntax)

# Install documentation as well:
add_subdirectory(doc)

# External projects...
include(setupExternalProjects)

# Setting up testing:
include(setupTesting)
# Add test folders:
include(ndifftests)
