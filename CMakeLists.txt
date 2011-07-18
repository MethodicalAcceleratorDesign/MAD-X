cmake_minimum_required(VERSION 2.6)

PROJECT(madX Fortran C CXX) 
#add name of languages used in the project as well

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
