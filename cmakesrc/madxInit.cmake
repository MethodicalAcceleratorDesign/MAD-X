###
#
# This file contains macros/functions for Mad-X, and
# sets some variables which are used by the system.
#
###


if ( MADX_FORCE_32 OR ${CMAKE_SIZEOF_VOID_P} EQUAL 4 )
    set(IS32BIT TRUE)
elseif (${CMAKE_SIZEOF_VOID_P} EQUAL 8)
    set(IS32BIT FALSE)
else()
    message(WARNING "Could not determine 32/64bit, assuming 32bit")
    set(IS32BIT TRUE)
endif()