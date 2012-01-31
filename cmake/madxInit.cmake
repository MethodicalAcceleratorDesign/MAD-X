###
#
# This file contains macros/functions for Mad-X, and
# sets some variables which are used by the system.
# 
# It also sets some initial stuff like versioning
###


if ( MADX_FORCE_32 OR ${CMAKE_SIZEOF_VOID_P} EQUAL 4 )
    set(IS32BIT TRUE)
elseif (${CMAKE_SIZEOF_VOID_P} EQUAL 8)
    set(IS32BIT FALSE)
else()
    message(WARNING "Could not determine 32/64bit, assuming 32bit")
    set(IS32BIT TRUE)
endif()


# project version
set( PROJECT_MAJOR_VERSION 5 )
set( PROJECT_MINOR_VERSION 00 )
set( PROJECT_PATCH_LEVEL 11 )

set( PROJECT_VERSION ${PROJECT_MAJOR_VERSION}.${PROJECT_MINOR_VERSION}.${PROJECT_PATCH_LEVEL} )


# Append _dev/-dev to binary/package name
if ( NOT PROJECT_PATCH_LEVEL EQUAL 00 )
    message(STATUS "Building a development version" )
    set (BINARY_POSTFIX "_dev")
    set (PKG_POSTFIX "-dev")
endif ( NOT PROJECT_PATCH_LEVEL EQUAL 00 )
