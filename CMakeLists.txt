# Suppress warnings on cygwin
# win32 is not anymore equal to 1 on cygwin
# with newer cmake. We still assume so though..
# Note:
# - In cygwin, the boolean CYGWIN is defined and 1
# - On Windows, WIN32 bool is defined, but 
# after cmake 2.8.4 it's normally zero when using
# cygwin.
set(CMAKE_LEGACY_CYGWIN_WIN32 0)

cmake_minimum_required(VERSION 2.6)

PROJECT(madX Fortran C CXX) 
#add name of languages used in the project as well

# Some policies regarding configure
# warnings we do not want..
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

# add source dir as a place for CMake modules (e.g. FindGSL.cmake)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmakesrc") 

# Setting project version:
include(setupVersion)

# Defining project options:
include(setupOptions)

# 
# we call a subscript that sets up our compiler specific flags...
# 
include(setupCompilerSpecifics)


# we call a subscript that will set up additional flags depending on 
# architecture chosen.
include(setupArchSpecifics)

#
# we call a subscript that sets up the source files:
# 
include(setupSources)

# Defining specific stuff for packaging:
# Should be called after sources are defined,
# before targets are defined.
include(setupPackaging)


# Setting up targets and linking between targets
# and external libraries:
include(setupTargets)

# 
# we call a subscript that sets up installation stuff:
# 
# Should be called after targets are defined
# 
include(setupInstallation)

# Setting up testing:
include(setupTesting)
