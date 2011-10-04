

###
#
# This file sets up specific flags for the PathScale compiler
#
###

if(CMAKE_Fortran_COMPILER_ID MATCHES "PathScale")
    set(CMAKE_Fortran_FLAGS_RELEASE " -funroll-loops -O3 ")
    if ( MADX_STATIC )
        set(CMAKE_Fortran_LINK_FLAGS   "${CMAKE_Fortran_LINK_FLAGS} -static ")
    endif()
else()
    message("You used setupPathScale.cmake but it had no effect")
endif()

