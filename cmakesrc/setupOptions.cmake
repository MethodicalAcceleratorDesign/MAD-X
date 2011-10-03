
# project options
if( NOT BUILD_SHARED_LIBS )
        option( BUILD_SHARED_LIBS "Turn on to build dynamic libraries"  OFF )
endif()

if (APPLE OR BUILD_SHARED_LIBS OR CYGWIN)
    option( MADX_STATIC "Turn on for static linking" OFF)
else()
    option ( MADX_STATIC "Turn on for static linking" ON)
endif()
# We need to specify what kind of library suffixes we search for in case
# for static linking:
if(MADX_STATIC)
    if(WIN32)
        set(CMAKE_FIND_LIBRARY_SUFFIXES .lib)
    elseif(APPLE)
        set(CMAKE_FIND_LIBRARY_SUFFIXES .a .lib)
    else()
        set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
    endif()
endif()

#Mad-X specific options (arch. specific options can be added in similar manner):
option( MADX_NTPSA "Build with NTPSA" ON)
option( MADX_FORCE_32 "Force 32bit build" OFF )
option( MADX_FEDORA_FIX "Fix for Fedora>11 for ifort compiler" OFF )

# double logic: first we set madx_online default on if sdds is found
# then if MADX_ONLINE is on without sdds found, we throw fatal error.

if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
#include our specific folders:
    if( MADX_FORCE_32 OR ${CMAKE_SIZEOF_VOID_P} EQUAL 4 )
        set(SDDS_SEARCH_DIRS  ${CMAKE_SOURCE_DIR}/lib32/)
    else()
        set(SDDS_SEARCH_DIRS  ${CMAKE_SOURCE_DIR}/lib64/)
    endif()
endif()
# normal search:
find_package(SDDS)

if(SDDS_FOUND)
    option( MADX_ONLINE "Build with Online model" ON)
else()
    option( MADX_ONLINE "Build with Online model" OFF)
endif()
if(MADX_ONLINE AND NOT SDDS_FOUND)
    message(FATAL_ERROR "SDDS is not found, required for the online model!")
endif()

# Default build type (defines different sets of flags)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
       "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel."
       FORCE)
endif()
