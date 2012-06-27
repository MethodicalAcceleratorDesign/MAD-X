
###
#
# This file sets up specific flags for the Intel compilers
#
###

if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -assume noold_unit_star -D_INTEL_IFORT_SET_RECL -fp-model precise -fpp")
    set(CMAKE_Fortran_FLAGS_RELEASE "-funroll-loops -O3")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O3 -g")
    set(CMAKE_Fortran_LINK_FLAGS   "${CMAKE_Fortran_LINK_FLAGS} -nofor_main")
    if ( MADX_STATIC )
        set(CMAKE_Fortran_LINK_FLAGS   "${CMAKE_Fortran_LINK_FLAGS} -static")
    endif ()
    if(IS32BIT)
        set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fp-model precise -m32")
    endif()
    set(INTEL_FLAGS "-Wcheck -Wp64 -diag-disable 2259,1572,981 -mp1 -fp-model strict")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x ${INTEL_FLAGS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${INTEL_FLAGS}")
    add_definitions(-D_ICC)
    add_definitions(-D_IFORT)
else()
    message("You used setupIntel.cmake but it had no effect")
endif()
