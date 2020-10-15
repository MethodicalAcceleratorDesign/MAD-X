# project options

option( MADX_STATIC "Turn on for static linking" OFF)
option( MADX_DEBUG "Turn on debug output" OFF)
option( USE_GC "Use Garbage Collector" ON)
option( MADX_NTPSA "Build with NTPSA" ON)
option( MADX_FORCE_32 "Force 32bit build" OFF)
option( MADX_LAPACK "Use system blas/lapack installation if one can be found" ON)

if (NOT (WIN32 OR CYGWIN))
 option( MADX_X11 "Turn on plotting using X11" ON )
endif()

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
    if (MADX_X11 AND (CMAKE_SYSTEM_NAME STREQUAL "Linux"))
       link_directories(${CMAKE_SOURCE_DIR}/lib${ARCH}/)    # for libX11.a
    endif()
endif()

#Mad-X specific options (arch. specific options can be added in similar manner):
if(APPLE)
 option( MADX_BUNDLE "Create bundle on OSX" OFF)
endif()

# double logic: first we set madx_online default on if sdds is found
# then if MADX_ONLINE is on without sdds found, we throw fatal error.

if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set(SDDS_SEARCH_DIRS  ${CMAKE_SOURCE_DIR}/lib${ARCH}/)
endif()
find_package(SDDS)

option( MADX_ONLINE "Build with Online model" ${SDDS_FOUND})
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

if(USE_GC)
   add_definitions("-D_USEGC")
endif()

if(MADX_DEBUG)
   add_definitions(-D_DEBUG -DDEBUG_ALL)
endif()

if(MADX_ONLINE)
   message(STATUS "Online model turned on")
   add_definitions(-D_ONLINE)
endif()
