###
#
# This file sets compiler specific flags
# and runs some special stuff needed for specific 
# compiler selections
# 
###

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${CMAKE_CURRENT_LIST_DIR}/compilers")

# These two compilers can be mixed at will...
include(setupGNU)
include(setupIntel)

# Only Fortran..
if(CMAKE_Fortran_COMPILER MATCHES "lf95")
   message( WARNING " This compiler is not supported for Mad-X")
   include(setupLahey)
elseif(CMAKE_Fortran_COMPILER MATCHES "nagfor")
   message( WARNING " This compiler is not supported for Mad-X")
   include(setupNAGFOR)
elseif(CMAKE_Fortran_COMPILER MATCHES "g77")
   message( WARNING " This compiler is not supported for Mad-X")
   message( "--- ifort is recommended fortran compiler ---")
   set(CMAKE_Fortran_FLAGS_RELEASE " -funroll-loops -fno-f2c -O3 ")
   set(CMAKE_Fortran_FLAGS_DEBUG   " -fno-f2c -O0 -g ")
   if ( MADX_STATIC )
      set(CMAKE_Fortran_LINK_FLAGS   "${CMAKE_Fortran_LINK_FLAGS} -static ")
   endif ()
elseif(CMAKE_Fortran_COMPILER MATCHES "g95")
   message( WARNING " This compiler is not supported for Mad-X")
   set(CMAKE_Fortran_FLAGS_RELEASE " -funroll-loops -fno-second-underscore -fshort-circuit -O2 ")
   set(CMAKE_Fortran_FLAGS_DEBUG   " -fno-second-underscore -O3 -g -Wall -pedantic -ggdb3")  
   add_definitions(-D_G95)
   if ( MADX_STATIC )
      set(CMAKE_Fortran_LINK_FLAGS   "${CMAKE_Fortran_LINK_FLAGS} -static ")
   endif ()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PathScale")
   message( WARNING " This compiler is not supported for Mad-X")
   include(setupPathScale)
endif()
#end fortran compiler stuff...


# General compile flags:
set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g ${CMAKE_Fortran_FLAGS_RELEASE}")
set(CMAKE_C_FLAGS_DEBUG   " ${CMAKE_C_FLAGS_DEBUG} -Wall -pedantic")
if(NOT (WIN32 AND CMAKE_C_COMPILER_ID STREQUAL "Intel"))
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -funroll-loops -std=c99")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -funroll-loops") #needed for c++ linking
set(CMAKE_CXX_FLAGS_DEBUG " ${CMAKE_CXX_FLAGS_DEBUG} -Wall")

add_definitions(-D_FULL)
add_definitions(-D_VERSION=${MADX_VERSION})
add_definitions(-D_VERSION_NUM=${VERSION_NUM})
add_definitions(-D_VERSION_DATE=${VERSION_DATE})
add_definitions(-D_VERSION_OSTYPE=${CMAKE_SYSTEM_NAME})

# C stuff:
# -- not needed for gnu/intel --
if(CMAKE_C_COMPILER_ID MATCHES "GNU" AND NOT CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
   execute_process(COMMAND ${C_COMPILER_NAME} -print-search-dirs
      OUTPUT_VARIABLE gccsearchdirs)
   string(REGEX REPLACE ".*libraries: =(.*)\n"  "\\1" gcclibs "${gccsearchdirs}")
   # need to do this many times because lf95 segfaults on lists with :
   string(REPLACE "/:/"  "/ -L/" gcclibs "${gcclibs}")
   # adding these to the linking process which is handled by a non-gnu fortran compiler in your case
   link_directories(${gcclibs}) 
endif()
# -- end of not needed for gnu/intel --
# end C stuff
