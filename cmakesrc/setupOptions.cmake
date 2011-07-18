
# project options
if( NOT BUILD_SHARED_LIBS )
        option( BUILD_SHARED_LIBS "Turn on to build dynamic libraries"  OFF )
endif()

if (APPLE OR BUILD_SHARED_LIBS)
    option( MADX_STATIC "Turn on for static linking" OFF)
else()
    option ( MADX_STATIC "Turn on for static linking" ON)
endif()

#Mad-X specific options (arch. specific options can be added in similar manner):
option( MADX_NTPSA "Build with NTPSA" ON)
option( MADX_FORCE_32 "Force 32bit build" OFF )
option( MADX_FEDORA_FIX "Fix for Fedora>11 for ifort compiler" OFF )

option( MADX_GOTOBLAS2 "Build with the GOTOBLAS2 libraries" OFF )
option( MADX_RICCARDO_FIX "Fix for Riccardo to find BLAS/LAPACK on his machine..." OFF )

if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    option( MADX_ONLINE "Build with Online model" ON)
else()
    option( MADX_ONLINE "Build with Online model" OFF)
endif()

# Default build type (defines different sets of flags)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
        "Choose the type of build, options are: Debug Release"
        FORCE)
endif()