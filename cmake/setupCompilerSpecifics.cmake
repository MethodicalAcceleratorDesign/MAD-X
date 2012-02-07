###
#
# This file sets compiler specific flags
# and runs some special stuff needed for specific 
# compiler selections
# 
###

if(CMAKE_VERSION VERSION_LESS "2.8.3")
  get_filename_component(CMAKE_CURRENT_LIST_DIR CMAKE_CURRENT_LIST_FILE PATH)
endif()
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${CMAKE_CURRENT_LIST_DIR}/compilers")

if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
   include(setupGNU)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    include(setupIntel)
elseif(CMAKE_Fortran_COMPILER MATCHES "lf95")
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
    if ( MADX_STATIC )
        set(CMAKE_Fortran_LINK_FLAGS   "${CMAKE_Fortran_LINK_FLAGS} -static ")
    endif ()

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PathScale")
    message( WARNING " This compiler is not supported for Mad-X")
    include(setupPathScale)

else()
    message( WARNING " This compiler is not supported for Mad-X")
    message("Fortran compiler full path: " ${CMAKE_Fortran_COMPILER})
    message("Fortran compiler: " ${Fortran_COMPILER_NAME})
    set(CMAKE_Fortran_FLAGS_RELEASE " -funroll-loops -fno-range-check -O2")
    if ( MADX_STATIC )
        set(CMAKE_Fortran_LINK_FLAGS   "${CMAKE_Fortran_LINK_FLAGS} -static ")
    endif ( MADX_STATIC )
endif()
#end fortran compiler stuff...


# General compile flags:
set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g ${CMAKE_Fortran_FLAGS_RELEASE}")
set(CMAKE_C_FLAGS_DEBUG   " ${CMAKE_C_FLAGS_DEBUG} -Wall -pedantic")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -funroll-loops -std=c99 -D_FULL ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -funroll-loops -std=c++98 -D_FULL ") #needed for c++ linking
set(CMAKE_CXX_FLAGS_DEBUG " ${CMAKE_CXX_FLAGS_DEBUG} -Wall")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -D_FULL ") 
if(MADX_DEBUG)
   add_definitions(-D_DEBUG -DDEBUG_ALL)
endif()

# Project version:
execute_process(COMMAND "date" "+%d-%m-%Y" OUTPUT_VARIABLE BUILD_DATE)
string(REGEX REPLACE "\n" "" BUILD_DATE ${BUILD_DATE})
add_definitions(-D_VERSION=${PROJECT_VERSION} -D_VERSION_DATE=${BUILD_DATE})

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

if(MADX_ONLINE)
    message(STATUS "Online model turned on")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_ONLINE ")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -D_ONLINE ")
endif()
