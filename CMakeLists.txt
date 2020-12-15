cmake_minimum_required(VERSION 2.8)

set(CMAKE_LEGACY_CYGWIN_WIN32 0) # Remove when CMake >= 2.8.4 is required

PROJECT(madX C CXX Fortran)

cmake_policy(VERSION 2.8...3.18)
if (NOT ${CMAKE_VERSION} VERSION_LESS 3.0)
    cmake_policy(SET CMP0042 NEW)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")
if(CMAKE_VERSION VERSION_LESS "2.8.3")
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/compilers")
endif()

include(madxInit)
include(setupOptions)

# We add libraries before setting up compilers because
# some of our options are incompatible with e.g. GC
if(USE_GC)
   set(CMAKE_POSITION_INDEPENDENT_CODE ON)
   add_subdirectory(libs/gc/gc-8.0.2)
endif()

include(setupCompilerSpecifics)
include(setupArchSpecifics)

add_subdirectory(libs/ptc)
add_subdirectory(libs/DISTlib)
add_subdirectory(src)
add_subdirectory(tools)
add_subdirectory(syntax)
add_subdirectory(doc)

include(setupTesting)
include(ndifftests)
