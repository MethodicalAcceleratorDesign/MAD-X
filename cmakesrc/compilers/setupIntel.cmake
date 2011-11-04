
###
#
# This file sets up specific flags for the Intel compilers
#
###

if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -assume noold_unit_star -D_INTEL_IFORT_SET_RECL -fp-model precise -fpp")
    set(CMAKE_Fortran_FLAGS_RELEASE "-funroll-loops -O3")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O3 -g")
    if ( MADX_STATIC )
        set(CMAKE_Fortran_LINK_FLAGS   "${CMAKE_Fortran_LINK_FLAGS} -static -nofor_main")
    endif ()
    if(IS32BIT)
        set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fp-model precise -m32")
    endif()
else()
    message("You used setupIntel.cmake but it had no effect")
endif()
