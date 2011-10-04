###
#
# This file sets up specific flags for the Lahey compiler
#
###


if(CMAKE_Fortran_COMPILER MATCHES "lf95")
    if ( MADX_FORCE_32 )
        message( WARNING " On a 64 bit system you need to use the toolchain-file (see README) to get anywhere with the 32bit compiler.")
    endif ( MADX_FORCE_32 )
    set(CMAKE_Fortran_FLAGS_RELEASE " --o2 --tp  ")
    set(CMAKE_SKIP_RPATH ON)
    set(CMAKE_Fortran_FLAGS_DEBUG   " --info --f95 --lst -V -g  --ap --trace --trap --verbose  --chk aesux ")
    set(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "") #suppress rdynamic which doesn't work for lf95...
    if ( MADX_STATIC )
        set(CMAKE_Fortran_LINK_FLAGS   "${CMAKE_Fortran_LINK_FLAGS} -static ")
    endif ()
    if(IS32BIT)
        set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wa,--32 --chk aesux ")
    else()
        set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} --chk aefs")
    endif()
else()
    message("You used setupLahey.cmake but it had no effect")
endif()
