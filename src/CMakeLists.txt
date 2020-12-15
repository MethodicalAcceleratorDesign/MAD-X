include(buildtypeCoverage)

file(GLOB src_files mad_*.c *.f90 *.F90 *.cpp)
file(GLOB main_file mad_main.c)

if(MADX_X11)
    file(GLOB gxx_remove gxx11ps.f90)
    set_property(SOURCE mad_gxx11c.c
                 APPEND PROPERTY COMPILE_DEFINITIONS _USE_X11)
else()
    file(GLOB gxx_remove gxx11.f90)
endif()

file(GLOB sdds_remove matchsdds.c)
file(GLOB mkthin_remove mad_mkthin.c)

if(MADX_LAPACK AND NOT MADX_FORCE_32)
   find_package(LAPACK)
endif()

if(LAPACK_FOUND AND BLAS_FOUND)
   message(STATUS "LAPACK uses ${LAPACK_LIBRARIES}")
   file(GLOB lapack_remove matchlib.f90 matchlib2.f90)
else()
   set_source_files_properties(matchlib2.f90 PROPERTIES COMPILE_FLAGS "-O0")
endif()

set(remove_list ${gxx_remove} ${sdds_remove} ${mkthin_remove} ${lapack_remove})
list(REMOVE_ITEM src_files ${remove_list})

file(GLOB headerfiles *.h)

add_library(madx ${src_files})

target_link_libraries(madx ptc)
target_link_libraries(madx DISTlib)

if(USE_GC)
   target_link_libraries(madx gc-lib)
endif()

# not sure if this is needed...
set_target_properties(madx PROPERTIES LINKER_LANGUAGE Fortran)

set_target_properties(madx     PROPERTIES VERSION ${MADX_VERSION})

if (${CMAKE_VERSION} VERSION_LESS "3.1")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
  set_target_properties(madx   PROPERTIES CXX_STANDARD 11)
endif()

if(MADX_BUNDLE)
   add_executable(madxbin MACOSX_BUNDLE  ${main_file})
else()
   add_executable(madxbin ${main_file})
endif()
set_target_properties(madxbin PROPERTIES LINKER_LANGUAGE Fortran VERSION ${MADX_VERSION})
set_target_properties(madxbin PROPERTIES OUTPUT_NAME "madx${BINARY_POSTFIX}")

target_link_libraries(madxbin madx)

if (MADX_X11)
    if (MADX_STATIC AND (CMAKE_SYSTEM_NAME STREQUAL "Linux"))
        target_link_libraries(madx X11)
    else()
        find_package(X11)
        if(X11_FOUND)
            include_directories(${X11_INCLUDE_DIR})
            target_link_libraries(madx ${X11_X11_LIB})
        endif()
    endif()
endif ()

if(LAPACK_FOUND AND BLAS_FOUND)
   target_link_libraries(madx ${LAPACK_LIBRARIES})
endif()

if(MADX_ONLINE)
   target_link_libraries(madx ${SDDS_LIBRARIES})
   target_link_libraries(madx z)
endif()

find_package(Threads)
find_package(STDCPP)
target_link_libraries(madx  ${CMAKE_THREAD_LIBS_INIT} ${STDCPP_LIBS})

if(USE_GC)
   include_directories(${CMAKE_SOURCE_DIR}/libs/gc/gc-8.0.2/include)
endif()

include(setupPackaging)
include(setupInstallation)
