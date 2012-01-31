
###
#
# This file sets up specific flags for the GNU compiler
#
###

if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
# General:
    set(CMAKE_Fortran_FLAGS " -fno-range-check -fno-f2c -cpp ") # remove -g -O2 from main list

    add_definitions(-D_GFORTRAN)
    # Release flags:
    # ON APPLE machines and on 32bit Linux systems, -O2 seems to be the highest optimization level possible
    # for file l_complex_taylor.f90
    if(APPLE OR ${CMAKE_SIZEOF_VOID_P} EQUAL 4)
        set(CMAKE_Fortran_FLAGS_RELEASE " -funroll-loops -O2 ")
    else()
        set(CMAKE_Fortran_FLAGS_RELEASE " -funroll-loops -O4 ")
    endif()

    # Additional option dependent flags:
    if ( MADX_STATIC )
        set(CMAKE_Fortran_LINK_FLAGS   "${CMAKE_Fortran_LINK_FLAGS} -static ")
    endif ()
    if(IS32BIT)
        set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -m32")
    endif()
else()
    message("You used setupGNU.cmake but it had no effect")
endif()
