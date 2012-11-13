
# project options

option( MADX_STATIC "Turn on for static linking" OFF)
option( MADX_DEBUG "Turn on debug output" OFF)

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
    if(BUILD_SHARED_LIBS)
       message(FATAL_ERROR "Cannot build shared libs with MADX_STATIC on")
    endif()
endif()

#Mad-X specific options (arch. specific options can be added in similar manner):
option( MADX_NTPSA "Build with NTPSA" ON)
option( MADX_FORCE_32 "Force 32bit build" OFF )
if(APPLE)
 option( MADX_BUNDLE "Create bundle on OSX" OFF)
endif()

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

if(MADX_NTPSA)
    add_definitions("-D_NTPSA")
endif()
