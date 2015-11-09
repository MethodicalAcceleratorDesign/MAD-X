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
file(READ ${CMAKE_SOURCE_DIR}/VERSION VERSION_INFO)
string(REGEX MATCH "VERSION *= *[0-9]+.[0-9]+.[0-9]+" MADX_VERSION  ${VERSION_INFO})
string(REGEX REPLACE "VERSION *= *" "" MADX_VERSION ${MADX_VERSION})
# split version in major/minor/patch_level:
string(REGEX REPLACE "^([0-9])+.[0-9]+.[0-9]+" "\\1" MADX_MAJOR_VERSION ${MADX_VERSION})
string(REGEX REPLACE "^[0-9]+.([0-9]+).[0-9]+" "\\1" MADX_MINOR_VERSION ${MADX_VERSION})
string(REGEX REPLACE "^[0-9]+.[0-9]+.([0-9]+)" "\\1" MADX_PATCH_LEVEL ${MADX_VERSION})
#VERSION_NUM:
string(REGEX MATCH "VERSION_NUM += *[^\n]+" VERSION_NUM  ${VERSION_INFO})
string(REGEX REPLACE "VERSION_NUM += *" "" VERSION_NUM  ${VERSION_NUM})
#VERSION_DATE:
string(REGEX MATCH "VERSION_DATE += *[^\n]+" VERSION_DATE  ${VERSION_INFO})
string(REGEX REPLACE "VERSION_DATE += *" "" VERSION_DATE  ${VERSION_DATE})

message(STATUS "Mad-X version: ${MADX_VERSION}")
message(STATUS "Version num: ${VERSION_NUM}")
message(STATUS "Version date: ${VERSION_DATE}")

# Append _dev/-dev to binary/package name
if(NOT MADX_PATCH_LEVEL EQUAL 00)
    message(STATUS "Building a development version")
    set (BINARY_POSTFIX "_dev")
    set (PKG_POSTFIX "-dev")
endif()

# add 32 to the name for 32bit binaries..
if(IS32BIT)
   set(BINARY_POSTFIX "${BINARY_POSTFIX}32")
endif()

# Location of fortran modules:
set(CMAKE_Fortran_MODULE_DIRECTORY
    ${PROJECT_BINARY_DIR}/include/fortran/madX CACHE PATH "Single Directory for all fortran modules."
)

