# adding library:
add_library(madx ${srcfiles})
# not sure if this is needed...
set_target_properties(madx PROPERTIES LINKER_LANGUAGE Fortran) 

# adding an executable:
add_executable(madxbin MACOSX_BUNDLE  mad_main.c)
set_target_properties(madxbin PROPERTIES LINKER_LANGUAGE Fortran)

# we want to add _dev to the binary in case this is a dev version for the automatic packaging...
set_target_properties(madxbin PROPERTIES OUTPUT_NAME "madx${BINARY_POSTFIX}")
get_target_property(binaryname madxbin OUTPUT_NAME)

# we need to link executable to our own library:
target_link_libraries(madxbin madx)

# I turn off search for libraries in case you are on Linux,
# to make sure we make use of the lib/lib64 folders
if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux")
    #necessary to search for X11 for OSX instead of directly including it
    find_package(X11)
    if(X11_FOUND)
        include_directories(${X11_INCLUDE_DIR})
        target_link_libraries(madx ${X11_X11_LIB})
    endif()
else()
    target_link_libraries(madx X11)
endif()

if(LAPACK_FOUND AND BLAS_FOUND)
    target_link_libraries(madx ${LAPACK_LIBRARIES})
endif()

# Online libraries:
if(MADX_ONLINE)
    target_link_libraries(madx ${SDDS_LIBRARIES}) 
endif()
target_link_libraries(madx z)
# new additions from Frank:
# On Debian Stable, stdc++ is not linked by default...
target_link_libraries(madx pthread c stdc++ gcc_eh)
