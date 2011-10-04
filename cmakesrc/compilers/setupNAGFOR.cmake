###
#
# This file sets up specific flags for the NAGFOR compiler
#
###

if(CMAKE_Fortran_COMPILER MATCHES "nagfor")
    message( WARNING " Make sure you use the same gcc as nagfor is compiled with, or linking WILL fail.")
    set(CMAKE_SKIP_RPATH ON)
    set(CMAKE_Fortran_FLAGS_RELEASE " -gline -maxcontin=100 -ieee=full -D_NAG ")
    set(CMAKE_Fortran_FLAGS_DEBUG   " -gline -maxcontin=100 -ieee=full -D_NAG -C=all -nan ")
    set(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "") #suppress rdynamic which isn't recognized by nagfor...
    if ( MADX_STATIC )
        set(CMAKE_Fortran_LINK_FLAGS   "${CMAKE_Fortran_LINK_FLAGS} -Bstatic ")
    endif ()
    if(IS32BIT)
        set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}  -Wc,-m32 -Wl,-m32 -abi=32")
    endif()
else()
    message("You used setupNAGFOR.cmake but it had no effect")
endif()
